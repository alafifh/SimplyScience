from AI_Func import * 

query = "Brain Cancer Research 2025"
pubmed_text = fetch_pubmed_abstracts(query)

chat = client.chats.create(model="models/gemini-flash-lite-latest")

extract_claims(pubmed_text, query, chat)

facts = get_facts()
for f in facts:
    print(f"{f['category']}: {f['text']}")
    if facts:
        first_id = facts[0]["id"]
        sources = get_sources(first_id)
        print("Sources:", sources)
