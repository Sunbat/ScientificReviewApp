import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import asyncio
import csv
from fpdf import FPDF
from SciRevApp import run_articles

articles = []  # Global list to store fetched articles

def on_submit():
    query = query_entry.get().strip()
    max_results_str = max_entry.get().strip()

    if not query:
        messagebox.showwarning("Input Error", "Please enter a search query.")
        return

    try:
        max_results = int(max_results_str) if max_results_str else 10
    except ValueError:
        messagebox.showwarning("Input Error", "Max results must be a number.")
        return

    submit_button.config(state=tk.DISABLED)
    progress.start()
    status_label.config(text="Fetching articles...")

    def run_async():
        global articles
        try:
            articles.clear()
            fetched = asyncio.run(run_articles(query, max_results))
            articles.extend(fetched)
            update_table()
            status_label.config(text=f"✅ Fetched {len(articles)} articles.")
        except Exception as e:
            messagebox.showerror("Error", str(e))
            status_label.config(text="❌ Failed to fetch articles.")
        finally:
            submit_button.config(state=tk.NORMAL)
            progress.stop()

    root.after(100, run_async)

def update_table():
    for row in table.get_children():
        table.delete(row)
    for art in articles:
        table.insert("", "end", values=(art['id'], art['year'], art['citations'], art['score']))

def export_csv():
    if not articles:
        messagebox.showwarning("No data", "Run a search first.")
        return
    path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
    if path:
        try:
            with open(path, "w", newline="", encoding='utf-8') as f:
                writer = csv.DictWriter(f, fieldnames=["id", "year", "citations", "score"])
                writer.writeheader()
                writer.writerows(articles)
            messagebox.showinfo("Export Complete", f"CSV saved to:\n{path}")
        except Exception as e:
            messagebox.showerror("Export Failed", str(e))

def export_pdf():
    if not articles:
        messagebox.showwarning("No data", "Run a search first.")
        return
    path = filedialog.asksaveasfilename(defaultextension=".pdf", filetypes=[("PDF files", "*.pdf")])
    if path:
        try:
            pdf = FPDF()
            pdf.add_page()
            pdf.set_font("Arial", size=12)
            pdf.cell(200, 10, txt="PubMed Article Relevance Scores", ln=True, align="C")
            pdf.ln(10)
            for art in articles:
                line = f"ID: {art['id']}, Year: {art['year']}, Citations: {art['citations']}, Score: {art['score']}"
                pdf.multi_cell(0, 10, txt=line)
            pdf.output(path)
            messagebox.showinfo("Export Complete", f"PDF saved to:\n{path}")
        except Exception as e:
            messagebox.showerror("Export Failed", str(e))

# ----------------- GUI SETUP -----------------

root = tk.Tk()
root.title("🔬 PubMed Scientific Review Tool")
root.geometry("800x600")

# --- Input Section ---
input_frame = tk.Frame(root)
input_frame.pack(pady=10)

tk.Label(input_frame, text="Search Query:").grid(row=0, column=0, padx=5, sticky="e")
query_entry = tk.Entry(input_frame, width=40)
query_entry.grid(row=0, column=1, padx=5)

tk.Label(input_frame, text="Max Results:").grid(row=0, column=2, padx=5, sticky="e")
max_entry = tk.Entry(input_frame, width=10)
max_entry.grid(row=0, column=3, padx=5)

submit_button = tk.Button(input_frame, text="Search PubMed", command=on_submit)
submit_button.grid(row=0, column=4, padx=10)

# --- Progress Bar & Status ---
progress = ttk.Progressbar(root, orient="horizontal", length=300, mode="indeterminate")
progress.pack(pady=5)
status_label = tk.Label(root, text="")
status_label.pack()

# --- Results Table ---
table_frame = tk.Frame(root)
table_frame.pack(expand=True, fill="both", padx=10, pady=10)

table = ttk.Treeview(table_frame, columns=("ID", "Year", "Citations", "Score"), show="headings")
for col in ("ID", "Year", "Citations", "Score"):
    table.heading(col, text=col)
    table.column(col, anchor="center", width=100 if col != "ID" else 200)
table.pack(side="left", fill="both", expand=True)

scrollbar = ttk.Scrollbar(table_frame, orient="vertical", command=table.yview)
scrollbar.pack(side="right", fill="y")
table.configure(yscrollcommand=scrollbar.set)

# --- Export Buttons ---
export_frame = tk.Frame(root)
export_frame.pack(pady=10)

tk.Button(export_frame, text="Export to CSV", command=export_csv).grid(row=0, column=0, padx=10)
tk.Button(export_frame, text="Export to PDF", command=export_pdf).grid(row=0, column=1, padx=10)

# --- Start GUI Loop ---
root.mainloop()
