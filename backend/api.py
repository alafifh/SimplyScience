from flask import Flask, request, jsonify
from flask_cors import CORS
from app import fetch_pubmed_abstracts, extract_claims, get_facts, get_sources, CLAIM_DB
from AI import query, max_results, client, chat

app = Flask(__name__)
CORS(app, origins=["https://alafifh.github.io"])

@app.route('/search', methods=['GET'])
def search():
    pubmed_text = fetch_pubmed_abstracts(query, max_results=max_results)

    chat = client.chats.create(model="models/gemini-flash-lite-latest")

    CLAIM_DB.clear()

    extract_claims(pubmed_text, query, chat)
    facts = get_facts()

    return jsonify(ok=True, facts=facts)

'''
@app.get("/health")
def health():
    return jsonify(ok=True)

@app.post("/search")
def search():
    data = request.get_json(silent=True) or {}
    query = (data.get("query") or "").strip()
    max_results = int(data.get("max_results") or 5)

    if not query:
        return jsonify(ok=False, error="No query provided"), 400

    pubmed_text = fetch_pubmed_abstracts(query, max_results=max_results)

    chat = client.chats.create(model="models/gemini-flash-lite-latest")

    CLAIM_DB.clear()

    extract_claims(pubmed_text, query, chat)
    facts = get_facts()

    return jsonify(ok=True, facts=facts)
'''