# Taxonomic hierarchy

Stanza version

* http://togostanza.org/stanza/lineage_information/help

Code

* https://github.com/togostanza/togostanza/blob/master/lineage_information_stanza/stanza.rb

Get text version of this API results

* curl -H "Accept: text/plain" "http://togogenome.org/sparqlist/api/lineage?tax=9606"


## Parameters

* `tax` taxonomy ID
  * default: 9606
  * example: 9606, 31033

## Endpoint

http://togogenome.org/sparql

## `result` SPARQL query

```sparql
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX taxo: <http://ddbj.nig.ac.jp/ontologies/taxonomy#>
PREFIX taxid: <http://identifiers.org/taxonomy/>
SELECT
  (REPLACE(STR(?tax), "http://identifiers.org/taxonomy/", "") AS ?tax_id)
  ?tax_label
  ?step
FROM <http://togogenome.org/graph/taxonomy>
WHERE {
  ?search_tax rdfs:label ?o FILTER (?search_tax = taxid:{{tax}}) .
  ?search_tax rdfs:subClassOf ?tax OPTION (transitive, t_direction 1, t_min(0), t_step("step_no") as ?step) .
  ?tax rdfs:label ?tax_label .
  FILTER (?tax != taxid:1)
}
ORDER BY DESC(?step)
```

## `output` JSON transformation

```javascript
({
  json({result}) {
    return result.results.bindings.map((row) => {
      const taxid = row.tax_id.value;
      const label = row.tax_label.value;
      return [taxid, label];
    });
  },
  text({result}) {
    return result.results.bindings.map((row) => {
      const taxid = row.tax_id.value;
      const label = row.tax_label.value;
      return label + " (" + taxid + ")"
    }).join(" / ");
  }
})
```
