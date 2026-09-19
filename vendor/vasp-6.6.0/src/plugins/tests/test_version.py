import re
from importlib.metadata import version
from pathlib import Path


def test_version():
    version_number = re.compile(r"(?P<part>major|minor|patch)\s*=\s*(?P<number>\d+)")
    with open(Path.cwd() / "../version.F", "r") as version_file:
        version_parts = {
            search_result["part"]: int(search_result["number"])
            for line in version_file
            if (search_result := version_number.search(line))
        }
    expected = (
        f"{version_parts['major']}.{version_parts['minor']}.{version_parts['patch']}"
    )
    assert version("vasp") == expected
