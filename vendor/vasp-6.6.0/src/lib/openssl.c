/*
 * OpenSSL wrapper functions
 *
 * To generate a private key you can use
 *    openssl genpkey -algorithm RSA -out private_key.pem -pkeyopt rsa_keygen_bits:2048
 * keep this file secret and do not commit it to version control!
 *
 * To extract the public key from the private key you can use
 *    openssl pkey -pubout -in private_key.pem -out public_key.pem
 * the public key can be shared and is used to verify signatures.
 */
#include <pwd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <openssl/opensslv.h>
#include <openssl/pem.h>
#include <openssl/rsa.h>

#define MIN_OPENSSL_VERSION 0x10100000L // OpenSSL 1.1.0
#if OPENSSL_VERSION_NUMBER < MIN_OPENSSL_VERSION
    #error "OpenSSL version is too old. Minimum required version is 1.1.0 (0x10100000L). Please update your OpenSSL header files and compiler environment."
#endif

const char* VASP_PUBLIC_KEY =
    "-----BEGIN PUBLIC KEY-----\n"
    "MIIBIjANBgkqhkiG9w0BAQEFAAOCAQ8AMIIBCgKCAQEAuR4GKVuFOFEQIEpsZIuO\n"
    "cDrkLYWDyEzdpC7qCj0HHZysFcajCfi5B7zWSKS51koX6dOG0mhIMVhcY36/L9/n\n"
    "BaAASJ8hqWeH5bPJvniUBY7fs+qxTzxWOf5tXL73iOpdJNbe7MWTgmyIX0n19L+y\n"
    "AXBhKBMXrAX12Fl3JTznHZCkBcVh9Ad9+0fm9gULkBW9xUYEfxfqDbYp+1EyQ9s2\n"
    "M9i9dxFEaKtwVyKJtr6ImZr7RqWojNapitb6l6GPddGKHd1gSrEKtvoipIWy2+km\n"
    "PdF1lmiBsH9MfE7LQuSSHsH9I753NW4aI22dRx4c2A9sSeRUIXy8+Cy6UI5B8neS\n"
    "EwIDAQAB\n"
    "-----END PUBLIC KEY-----\n";


struct string_with_length
{
    char* string;
    int length;
};

/* Decode string from base64 to byte */
void base64_decode(
    const char* encoded_string, unsigned char** decoded_string, size_t* length
)
{
    BIO* bio = BIO_new_mem_buf(encoded_string, -1);
    BIO* b64 = BIO_new(BIO_f_base64());
    bio = BIO_push(b64, bio);

    int decode_length = strlen(encoded_string);
    *decoded_string = (unsigned char*)malloc(decode_length);
    if (*decoded_string)
    {
        *length = BIO_read(bio, *decoded_string, decode_length);
    }
    else
    {
        *length = 0;
    }
    BIO_free_all(bio);
}

/* Encode string from byte to base64 */
void base64_encode(
    const unsigned char* decoded_string, size_t length, char** encoded_string
)
{
    BUF_MEM *buffer_ptr;
    BIO* b64 = BIO_new(BIO_f_base64());
    BIO* bio = BIO_new(BIO_s_mem());
    bio = BIO_push(b64, bio);

    BIO_set_flags(bio, BIO_FLAGS_BASE64_NO_NL); // No newlines
    BIO_write(bio, decoded_string, length);
    BIO_flush(bio);
    BIO_get_mem_ptr(bio, &buffer_ptr);
    BIO_set_close(bio, BIO_NOCLOSE);
    BIO_free_all(bio);

    strncpy(*encoded_string, buffer_ptr->data, buffer_ptr->length);
    (*encoded_string)[buffer_ptr->length] = '\0'; // Null-terminate the string
    BUF_MEM_free(buffer_ptr);
}

/* Convert public key from PEM format to internal OpenSSL format */
EVP_PKEY* create_evp_from_public_key(const char* public_key)
{
    BIO* bio = BIO_new_mem_buf((void*)public_key, -1);
    EVP_PKEY* evp_pkey = PEM_read_bio_PUBKEY(bio, NULL, NULL, NULL);
    BIO_free_all(bio);
    return evp_pkey;
}


/* Verify the signature of a message using the public key */
int verify_signature(
    const struct string_with_length* message,
    const struct string_with_length* signature_encoded
)
{
    unsigned char* signature_decoded = NULL;
    size_t signature_length = 0;
    int error = -1;  // Error decoding signature
    base64_decode(signature_encoded->string, &signature_decoded, &signature_length);
    if (!signature_decoded) goto finalize;

    error = -2;  // Error creating EVP_PKEY
    EVP_PKEY* evp_pkey = create_evp_from_public_key(VASP_PUBLIC_KEY);
    if (!evp_pkey) goto finalize;

    error = -3;  // Error creating EVP_MD_CTX
    EVP_MD_CTX* ctx = EVP_MD_CTX_new();
    if (!ctx) goto finalize;

    error = -4;  // Error initializing DigestVerify
    if (EVP_DigestVerifyInit(ctx, NULL, EVP_sha256(), NULL, evp_pkey) <= 0) goto finalize;

    error = -5;  // Error during DigestVerifyUpdate
    if (EVP_DigestVerifyUpdate(ctx, message->string, message->length) <= 0) goto finalize;

    int verify_status = EVP_DigestVerifyFinal(ctx, signature_decoded, signature_length);
    if (verify_status == 1)
    {
        error = 0;
    }
    else if (verify_status == 0)
    {
        error = 1; // Signature is invalid
    }
    else
    {
        error = -6; // Some other error occurred during verification
    }

finalize:
    if (signature_decoded) free(signature_decoded);
    if (ctx) EVP_MD_CTX_free(ctx);
    if (evp_pkey) EVP_PKEY_free(evp_pkey);

    return error;
}


/* Get user name using getpwuid() */
int get_user_name(struct string_with_length* user_name)
{
    struct passwd* pw = getpwuid(getuid());
    if (!pw || !pw->pw_name)
    {
        user_name->length = 0;
        return -7; // Error: could not get user name
    }
    user_name->length = strlen(pw->pw_name);
    strncpy(user_name->string, pw->pw_name, user_name->length);
    return 0;
}

/* Get user id from syscall */
long get_user_id()
{
    return syscall(SYS_getuid);
}

/* Hash user name and user id with SHA-256 */
void hash_user_data(const struct string_with_length* user_data, struct string_with_length* hash)
{
    unsigned char digest[SHA256_DIGEST_LENGTH];
    SHA256((unsigned char*)user_data->string, user_data->length, digest);
    base64_encode(digest, SHA256_DIGEST_LENGTH, &hash->string);
    hash->length = strlen(hash->string);
}
