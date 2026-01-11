from app import * 

query = "effects of intermittent fasting on metabolic health"
pubmed_text = fetch_pubmed_abstracts(query)

print ()
print ()
extract_claims(pubmed_text, query)
print ()
print ()
facts = get_facts()
for f in facts:
    print(f"{f['category']}: {f['text']}")
    if facts:
        first_id = facts[0]["id"]
        sources = get_sources(first_id)
        print(sources)
print ()
print ()