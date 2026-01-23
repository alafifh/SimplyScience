from app import fetch_pubmed_abstracts, extract_claims, get_facts, get_sources, chat

query = "gauss's law"
pubmed_text = fetch_pubmed_abstracts(query)

print ()
print ()
extract_claims(pubmed_text, query, chat)
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