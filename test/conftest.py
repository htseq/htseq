import os
import sysconfig
import pytest


def get_data_folder():
    cwd = os.getcwd().rstrip("/")
    if cwd.endswith("example_data"):
        return cwd + "/"
    else:
        return cwd + "/example_data/"


def get_docs_folder():
    cwd = os.getcwd().rstrip("/")
    if cwd.endswith("doc"):
        return cwd + "/"
    elif cwd.endswith("example_data"):
        return cwd[: -len("example_data")] + "doc/"
    else:
        return cwd + "/doc/"


def get_scripts_folder():
    return sysconfig.get_path("scripts") + "/"


# Same as fixtures
@pytest.fixture(scope="module")
def data_folder():
    return get_data_folder()


@pytest.fixture(scope="module")
def docs_folder():
    return get_docs_folder()


@pytest.fixture(scope="module")
def scripts_folder():
    return get_scripts_folder()


def _clean_Python_warnings(text):
    """Clean Python warnings from a stderr call.

    Parameters:
        text (str): Text to clean.
    Returns:
        str: Cleaned text.
    """
    new_lines = []
    for line in text.split("\n"):
        if line.startswith("<") and "Warning:" in line:
            continue
        new_lines.append(line)
    return "\n".join(new_lines)
