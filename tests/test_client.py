"""uniprot_client: pagination link parsing (no network)."""
from pipeline.uniprot_client import _NEXT_LINK


def test_next_link_regex():
    hdr = '<https://rest.uniprot.org/uniprotkb/search?cursor=abc&query=x&size=500>; rel="next"'
    assert _NEXT_LINK.search(hdr).group(1).startswith("https://rest.uniprot.org/uniprotkb/search?cursor=abc")
    assert _NEXT_LINK.search("") is None
