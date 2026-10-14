from pathlib import Path
import re
from urllib.parse import unquote, urlsplit

import pytest


ROOT = Path(__file__).resolve().parents[1]
DOCUMENTS = [ROOT / "README.md", ROOT / "technical_note.md", *sorted((ROOT / "reports").glob("*/*.md"))]


@pytest.mark.parametrize("document", DOCUMENTS, ids=lambda path: str(path.relative_to(ROOT)))
def test_relative_documentation_links_resolve(document):
    targets = re.findall(r"!?\[[^\]]*\]\(([^)]+)\)", document.read_text())
    for target in targets:
        target = target.split(' "', 1)[0].strip("<>")
        parsed = urlsplit(target)
        if parsed.scheme or parsed.netloc or not parsed.path:
            continue
        assert (document.parent / unquote(parsed.path)).exists(), f"{document.relative_to(ROOT)}: {target}"
