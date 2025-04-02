from Bio import Entrez
from datetime import datetime
import re

# Set your email here for Entrez API usage
email = input("Enter your email: ")
Entrez.email = email

def search_pubmed(query, max_results=10):
    """Search PubMed for articles matching the query."""
    handle = Entrez.esearch(db='pubmed', 
                             sort='relevance', 
                             retmax=max_results,
                             term=query)
    results = Entrez.read(handle)
    handle.close()
    return results['IdList']

def get_article_details(article_id):
    """Retrieve article details from PubMed."""
    handle = Entrez.efetch(db='pubmed', id=article_id, rettype='medline', retmode='text')
    article_details = handle.read()
    handle.close()
    return article_details

def extract_publication_date(article_details):
    """Extract publication date from article details."""
    match = re.search(r'DP\s+-\s+(\d{4})', article_details)
    if match:
        return match.group(1)
    return "0000"

def extract_citation_count(article_details):
    """Extract citation count from article details (if available)."""
    # Citation count might not be in the retrieved data; this is an example regex.
    # Adjust based on actual data format.
    match = re.search(r'VI\s+-\s+(\d+)', article_details)
    if match:
        return int(match.group(1))
    return 0
# update relevance scoring algorithm asap
# elaborate on calculation algo
def calculate_relevance_score(publication_year, citation_count):
    """Calculate the relevance score based on publication year and citation count."""
    current_year = datetime.now().year
    years_since_publication = current_year - int(publication_year)
    relevance_score = citation_count * (1 / (years_since_publication + 1))
    return relevance_score

def main():
    """Main function to execute the PubMed search and relevance scoring."""
    query = input("Enter your PubMed search query: ")
    max_results = int(input("Enter the maximum number of results to retrieve: "))
    
    article_ids = search_pubmed(query, max_results)
    print(f"Retrieving information for {len(article_ids)} articles...\n")
    
    for article_id in article_ids:
        article_details = get_article_details(article_id)
        publication_year = extract_publication_date(article_details)
        citation_count = extract_citation_count(article_details)
        relevance_score = calculate_relevance_score(publication_year, citation_count)
        
        # Print article information
        print(f"PubMed ID: {article_id}")
        print(f"Publication Year: {publication_year}")
        print(f"Citation Count: {citation_count}")
        print(f"Relevance Score: {relevance_score}\n")

if __name__ == "__main__":
    main()

# Review relevance scoring algo
# UI?
# Give score info and reference
# Relevance scoring algo, improvement
    #