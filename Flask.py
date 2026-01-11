from flask import Flask, render_template, request, jsonify
import json

app = Flask(__name__)

@app.route('/')
def index():
    return render_template("searchtest.html")

@app.route('/search', methods=['POST'])
def search():
    query = request.json.get("query")
    if not query:
        return jsonify({"error": "No search term provided"}), 400

    # Fetch abstracts
    abstracts = fetch_pubmed_abstracts(query, max_results=5)

    # Create a Gemini chat
    chat = client.chats.create(model="models/gemini-flash-lite-latest")
    
    # Extract claims into CLAIM_DB
    extract_claims(abstracts, query, chat)

    # Return facts for frontend display
    facts = get_facts()
    return jsonify(facts)

@app.route('/sources/<claim_id>')
def sources(claim_id):
    return jsonify({"sources": get_sources(claim_id)})

if __name__ == "__main__":
    app.run(debug=True)
