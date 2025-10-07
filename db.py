import sqlite3
from flask import g
   
    
DATABASE = "asteroids.db"


# --- DB init ---
def init_db():
    with sqlite3.connect(DATABASE) as db, open("schema.sql", "r") as f:
        db.executescript(f.read())
        
        
# --- DB conn ---
def get_db():
    if "db" not in g:
        g.db = sqlite3.connect(DATABASE)
        g.db.row_factory = sqlite3.Row
    return g.db


def close_db(exception):
    db = g.pop("db", None)
    if db is not None:
        db.close()
