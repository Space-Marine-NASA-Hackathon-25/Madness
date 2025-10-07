from flask import Flask
from flask_session import Session
from db import init_db, close_db
from routes import bp as routes_bp


def create_app():
    app = Flask(__name__)
    app.register_blueprint(routes_bp)
    
    app.config["SESSION_PERMANENT"] = False
    app.config["SESSION_TYPE"] = "filesystem"
    Session(app)

    app.teardown_appcontext(close_db)
    
    with app.app_context():
        init_db()      
         
    return app


if __name__ == "__main__": 
    app = create_app()   
    app.run(debug=True, port=5000)