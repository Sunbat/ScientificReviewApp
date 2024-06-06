from Bio import Entrez
# from pymed import PubMed
from datetime import datetime

def search_pubmed(query, max_results=10):
    Entrez.email = 'loicrulin@gmail.com'  # Enter your email here
    handle = Entrez.esearch(db='pubmed', 
                             sort='relevance', 
                             retmax=max_results,
                             term=query)
    results = Entrez.read(handle)
    handle.close()
    return results['IdList']

def get_article_details(article_id):
    Entrez.email = 'your_email@example.com'  # Enter your email here
    handle = Entrez.efetch(db='pubmed', id=article_id, rettype='medline', retmode='text')
    article_details = handle.read()
    handle.close()
    return article_details

def calculate_relevance_score(article_date, citation_count):
    # Calculate the age of the article in years
    article_date = datetime.strptime(article_date, '%Y%m%d')
    years_since_publication = (datetime.now() - article_date).days / 365
    
    # Calculate relevance score
    relevance_score = citation_count * (1 / (years_since_publication + 1))
    return relevance_score

def main():
    query = input("Enter your PubMed search query: ")
    max_results = int(input("Enter the maximum number of results to retrieve: "))
    
    article_ids = search_pubmed(query, max_results)
    print(f"Retrieving information for {len(article_ids)} articles...\n")
    
    for article_id in article_ids:
        article_details = get_article_details(article_id)
        # Parse article details to get publication date and citation count
        article_date = ''  # Extract publication date
        citation_count = 0  # Extract citation count
        relevance_score = calculate_relevance_score(article_date, citation_count)
        
        # recheck git doc
        
        # Print article information
        print(f"PubMed ID: {article_id}")
        print(f"Publication Date: {article_date}")
        print(f"Citation Count: {citation_count}")
        print(f"Relevance Score: {relevance_score}\n")

if __name__ == "__main__":
    main()
