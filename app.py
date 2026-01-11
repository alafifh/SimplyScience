from flask import Flask, request, jsonify
from AI_Func import * 

app = Flask(__name__)

@app.route('/search', methods=['POST'])
def search():
    data = request.get_json()
    query = data.get('query')
    if not query:
        return jsonify({"error": "No query provided"}), 400

    pubmed_text = fetch_pubmed_abstracts(query)
    extract_claims(pubmed_text, query, chat)  # `chat` = your Gemini client chat

    facts = get_facts()  # Return all high/moderate evidence facts
    return jsonify(facts)