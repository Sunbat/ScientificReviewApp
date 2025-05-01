import asyncio
from SciRevApp import run_articles

def display_results(articles):
    print("\n📄 Top PubMed Articles by Relevance:\n")
    print(f"{'ID':<12} {'Year':<6} {'Citations':<10} {'Score':<8} Title")
    print("-" * 100)
    for art in articles:
        print(f"{art['id']:<12} {art['year']:<6} {art['citations']:<10} {art['score']:<8} {art['title']}")
    print("\n✅ Done.")

def main():
    print("🔬 PubMed Scientific Review (CLI Tool)")

    query = input("Enter your PubMed search query (e.g., cancer genomics): ").strip()
    if not query:
        print("❌ Query cannot be empty.")
        return

    try:
        max_results = int(input("Max number of results [default = 10]: ").strip() or "10")
    except ValueError:
        print("Invalid input. Using default = 10")
        max_results = 10

    print("\n🔎 Searching PubMed... This may take a few seconds...\n")
    articles = asyncio.run(run_articles(query, max_results))

    if not articles:
        print("⚠️ No articles found or an error occurred.")
    else:
        display_results(articles)

if __name__ == "__main__":
    main()
