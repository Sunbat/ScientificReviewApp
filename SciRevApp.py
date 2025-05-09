import re
from datetime import datetime
from Bio import Entrez
import aiohttp
import asyncio
import urllib.parse
import math

# --- Ask user for PubMed email once ---
Entrez.email = input("📧 Enter your email for PubMed (Entrez API): ").strip()
if not Entrez.email:
    raise ValueError("Email is required to use PubMed Entrez API.")

# --- PubMed Search ---
async def search_pubmed(query, max_results=10):
    handle = Entrez.esearch(db='pubmed', sort='relevance', retmax=max_results, term=query)
    results = Entrez.read(handle)
    handle.close()
    return results["IdList"]

# --- Fetch Article Details, Title, and DOI ---
async def fetch_article_detail(article_id):
    try:
        # Fetch full MEDLINE text
        handle = Entrez.efetch(db='pubmed', id=article_id, rettype='medline', retmode='text')
        details = handle.read()
        handle.close()

        # Try to extract DOI from raw text (fallback if summary fails)
        doi_match = re.search(r'(10\.\d{4,9}/[-._;()/:A-Z0-9]+)', details, re.IGNORECASE)
        regex_doi = doi_match.group(1) if doi_match else ""

        # Fetch structured metadata for title and DOI
        summary_handle = Entrez.esummary(db='pubmed', id=article_id)
        summary = Entrez.read(summary_handle)
        summary_handle.close()

        title = summary[0].get("Title", "")
        summary_doi = summary[0].get("DOI", "")

        doi = summary_doi or regex_doi

        return article_id, details, doi, title
    except Exception as e:
        print(f"❌ Error retrieving article {article_id}: {e}")
        return article_id, "", "", ""

# --- Year Extraction from MEDLINE ---
def extract_publication_year(details):
    match = re.search(r'DP\s+-\s+(\d{4})', details)
    return int(match.group(1)) if match else 2000

# --- Match query terms against MEDLINE text ---
def compute_term_match_score(query, details):
    query_terms = query.lower().split()
    hits = sum(1 for word in query_terms if word in details.lower())
    return min(hits / len(query_terms), 1.0)

# --- CrossRef fallback if DOI missing ---
async def find_doi_by_title(title):
    if not title:
        return ""
    query = urllib.parse.quote(title)
    url = f"https://api.crossref.org/works?query.title={query}&rows=1"
    try:
        async with aiohttp.ClientSession() as session:
            async with session.get(url) as resp:
                if resp.status == 200:
                    data = await resp.json()
                    items = data.get("message", {}).get("items", [])
                    if items and "DOI" in items[0]:
                        return items[0]["DOI"]
    except Exception as e:
        print(f"⚠️ CrossRef error for title '{title}': {e}")
    return ""

# --- OpenCitations query ---
async def extract_citation_count(doi):
    if not doi:
        return 0

    encoded_doi = urllib.parse.quote(doi)
    opencitations_url = f"https://opencitations.net/index/coci/api/v1/citations/{encoded_doi}"
    semantic_url = f"https://api.semanticscholar.org/graph/v1/paper/DOI:{encoded_doi}?fields=citationCount"

    try:
        async with aiohttp.ClientSession() as session:
            # Try OpenCitations first
            async with session.get(opencitations_url) as resp:
                if resp.status == 200:
                    data = await resp.json()
                    if isinstance(data, list) and len(data) > 0:
                        return len(data)

            # Fallback: try Semantic Scholar
            async with session.get(semantic_url) as resp:
                if resp.status == 200:
                    data = await resp.json()
                    return data.get("citationCount", 0)
    except Exception as e:
        print(f"⚠️ Citation fetch error for DOI {doi}: {e}")

    return 0

# --- Relevance Score ---
def calculate_relevance_score(year, citations, term_score):
    current_year = datetime.now().year
    years_since = current_year - year
    if years_since < 0:  # Handles potential future dates, though unlikely for publications
        years_since = 0
    log_citations = math.log10(citations + 1)
    normalization_factor_citations = math.log10(1000 + 1)
    if normalization_factor_citations == 0:  # Should not happen with 1000+1
        citation_score_component = 0.0
    else:
        citation_score_component = log_citations / normalization_factor_citations

    citation_score = min(citation_score_component, 1.0)  # Cap the score at 1.0
    w1, w2, w3 = 0.4, 0.3, 0.3
    recency_score = 1 / (1 + years_since)
    final_score = (w1 * citation_score) + \
                  (w2 * recency_score) + \
                  (w3 * term_score)

    return round(final_score, 4)

# --- Main Entry Point for GUI/CLI/Web ---
async def run_articles(query, max_results=10):
    article_ids = await search_pubmed(query, max_results)
    fetch_tasks = [fetch_article_detail(aid) for aid in article_ids]
    fetched = await asyncio.gather(*fetch_tasks)

    articles = []

    for aid, details, doi, title in fetched:
        year = extract_publication_year(details)
        term_score = compute_term_match_score(query, details)

        if not doi:
            doi = await find_doi_by_title(title)

        print(f"🔍 Article ID: {aid} | DOI: {doi or '❌ Not Found'}")

        citations = await extract_citation_count(doi)
        score = calculate_relevance_score(year, citations, term_score)

        articles.append({
            "id": aid,
            "year": year,
            "citations": citations,
            "score": score,
            "title": title
        })

    return articles
