from AI_Func import * 

query = "effects of intermittent fasting on metabolic health"
pubmed_text = fetch_pubmed_abstracts(query)

chat = client.chats.create(model="models/gemini-flash-lite-latest")
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
