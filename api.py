from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

from AI_Func import (
    client,
    fetch_pubmed_abstracts,
    extract_claims,
    get_facts,
    CLAIM_DB,
)

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["https://alafifh.github.io"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# -------- Request schema --------
class FeedRequest(BaseModel):
    query: str
    max_results: int = 5

@app.get("/health")
def health():
    return {"ok": True}

@app.post("/feed")
def feed(req: FeedRequest):
    """
    USER provides the query.
    """
    query = req.query.strip()
    if not query:
        raise HTTPException(status_code=400, detail="Query is required")

    # 1) Fetch PubMed abstracts
    pubmed_text = fetch_pubmed_abstracts(query, max_results=req.max_results)

    # 2) Create Gemini chat session
    chat = client.chats.create(
        model="models/gemini-flash-lite-latest"
    )

    # 3) Clear old claims (VERY important)
    CLAIM_DB.clear()

    # 4) Run AI extraction
    extract_claims(pubmed_text, query, chat)

    # 5) Return structured facts
    facts = get_facts()
    return {
        "ok": True,
        "query": query,
        "facts": facts
    }
