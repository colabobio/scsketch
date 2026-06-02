from scsketch._api import fetch_gene_description


class _Response:
    def __init__(self, payload, status_code=200):
        self._payload = payload
        self.status_code = status_code

    def raise_for_status(self):
        if self.status_code >= 400:
            raise AssertionError(f"unexpected status {self.status_code}")

    def json(self):
        return self._payload


def test_fetch_gene_description_queries_symbol(monkeypatch):
    fetch_gene_description.cache_clear()
    calls = []

    def fake_get(url, *, params=None, timeout=None):
        calls.append((url, params, timeout))
        return _Response(
            {
                "hits": [
                    {
                        "_id": "7157",
                        "symbol": "TP53",
                        "name": "tumor protein p53",
                        "summary": "Acts as a tumor suppressor.",
                    }
                ]
            }
        )

    monkeypatch.setattr("scsketch._api.requests.get", fake_get)

    result = fetch_gene_description("TP53")

    assert result["symbol"] == "TP53"
    assert result["summary"] == "Acts as a tumor suppressor."
    assert calls == [
        (
            "https://mygene.info/v3/query",
            {
                "q": "symbol:TP53",
                "fields": "symbol,name,summary,entrezgene,ensembl.gene,taxid",
                "species": "human",
                "size": 1,
            },
            15,
        )
    ]


def test_fetch_gene_description_uses_cache(monkeypatch):
    fetch_gene_description.cache_clear()
    calls = []

    def fake_get(url, *, params=None, timeout=None):
        calls.append((url, params, timeout))
        return _Response({"hits": [{"symbol": "CDK2", "name": "CDK2"}]})

    monkeypatch.setattr("scsketch._api.requests.get", fake_get)

    assert fetch_gene_description("CDK2")["symbol"] == "CDK2"
    assert fetch_gene_description("CDK2")["symbol"] == "CDK2"
    assert len(calls) == 1


def test_fetch_gene_description_fetches_ensembl_annotation(monkeypatch):
    fetch_gene_description.cache_clear()
    calls = []

    def fake_get(url, *, params=None, timeout=None):
        calls.append((url, params, timeout))
        return _Response(
            {
                "_id": "ENSG00000141510",
                "symbol": "TP53",
                "name": "tumor protein p53",
                "summary": "Acts as a tumor suppressor.",
            }
        )

    monkeypatch.setattr("scsketch._api.requests.get", fake_get)

    result = fetch_gene_description("ENSG00000141510")

    assert result["symbol"] == "TP53"
    assert calls == [
        (
            "https://mygene.info/v3/gene/ENSG00000141510",
            {"fields": "symbol,name,summary,entrezgene,ensembl.gene,taxid"},
            15,
        )
    ]


def test_fetch_gene_description_returns_none_without_match(monkeypatch):
    fetch_gene_description.cache_clear()

    def fake_get(url, *, params=None, timeout=None):
        return _Response({"hits": []})

    monkeypatch.setattr("scsketch._api.requests.get", fake_get)

    assert fetch_gene_description("NOT_A_GENE") is None
