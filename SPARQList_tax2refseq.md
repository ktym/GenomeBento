# Taxonomy to RefSeq

Stanza version

* http://togostanza.org/stanza/genome_information/help

Code

* https://github.com/togostanza/togostanza/blob/master/genome_information_stanza/stanza.rb

Get FASTA sequences via TogoWS

* curl -H "Accept: text/plain" "http://togogenome.org/sparqlist/api/tax2refseq?tax=9606" | sh

## Parameters

* `tax` taxonomy ID
  * default: 9606
  * example: 9606, 31033

## Endpoint

http://togogenome.org/sparql

## `result` SPARQL query

```sparql
DEFINE sql:select-option "order"
PREFIX obo: <http://purl.obolibrary.org/obo/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX insdc: <http://ddbj.nig.ac.jp/ontologies/nucleotide/>
PREFIX idtax: <http://identifiers.org/taxonomy/>
PREFIX togo: <http://togogenome.org/stats/>
PREFIX ass: <http://ddbj.nig.ac.jp/ontologies/assembly/>
SELECT DISTINCT ?refseq_version ?taxid
FROM <http://togogenome.org/graph/refseq>
FROM <http://togogenome.org/graph/so>
FROM <http://togogenome.org/graph/stats>
FROM <http://togogenome.org/graph/assembly_report>
WHERE
{
  VALUES ?taxid { "{{tax}}" }
  ?seq obo:RO_0002162 idtax:{{tax}} .
  ?refseq_link insdc:sequence ?seq ;
    a insdc:Entry ;
    insdc:sequence_version ?refseq_version .
}
ORDER BY ?refseq_version
```

## `output` JSON transformation

```javascript
({
  json({result}) {
    return result.results.bindings.map((row) => {
      return row.refseq_version.value;
    });
  },
  text({result}) {
    return result.results.bindings.map((row) => {
      const refseq = row.refseq_version.value;
      const taxid = row.taxid.value;
      return "curl http://togows.org/entry/nucleotide/" + refseq + ".fasta >> " + taxid + ".fasta";
    }).join("\n");
  }
})
```
