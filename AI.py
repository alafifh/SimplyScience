from AI_Func import * 

if __name__ == "__main__":
    query = "Brain Cancer Research 2025"
    pubmed_text = fetch_pubmed_abstracts(query)

    # Create ONE Gemini chat session (free-tier)
    chat = client.chats.create(
        model="models/gemini-flash-lite-latest"
    )

    # Extract structured claims
    extract_claims(pubmed_text, query, chat)

    # Display facts
    print(query)
    facts = get_facts()
    for f in facts:
        print(f"{f['category']}: {f['text']} (Evidence: {f['evidence']})\n")

    # Example: show sources for first claim
    if facts:
        first_id = facts[0]["id"]
        sources = get_sources(first_id)
        print("=== Sources for this claim ===")
        print(sources)
