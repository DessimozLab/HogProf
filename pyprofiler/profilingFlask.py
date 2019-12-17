#run a profiler as a flask service

#send back info on returned hogs, linkage, tax spread

#

import profiler
from utils import config_utils
from flask import jsonify
import argparse
from flask import Flask, make_response
app = Flask(__name__)
from flask import Blueprint

bp = Blueprint('api', __name__)

from app.api import users, errors, tokens


@app.route('/knn/<str:protid>')
def return_knn(protid):
    #return csv data of the output matrix

    response = make_response(csv)
    cd = 'attachment; filename=mycsv.csv'
    response.headers['Content-Disposition'] = cd
    response.mimetype='text/csv'
    return response

@app.route('/profile/<str:protid>')
def return_profile():
    #return csv data of the output matrix
    response = make_response(csv)
    return response


@app.dist('/profile/<str:protids>')
def download_csv(protids):
    #return distance between two profiles
    csv = 'foo,bar,baz\nhai,bai,crai\n'
    response = make_response(csv)
    cd = 'attachment; filename=mycsv.csv'
    response.headers['Content-Disposition'] = cd
    response.mimetype='text/csv'
    return response


@bp.route('/users/<int:id>', methods=['GET'])
def get_user(id):
    pass

@bp.route('/profiler/<str:id>/knn', methods=['GET'])
def get_followers(id):
    pass

@bp.route('/profiler/<str:id>/profile', methods=['GET'])
def get_profile(id):
    pass


@bp.route('/users/<int:id>/followed', methods=['GET'])
def get_followed(id):
    pass

@bp.route('/users', methods=['POST'])
def get_distmat():
    data = request.get_json() or {}
    #send array result as json output
    
    response = jsonify(array2lists(distmat) )

    return response

@bp.route('/users/<int:id>', methods=['PUT'])
def update_user(id):
    pass

def bad_request():
    pass

def get_token():
    pass

def revoke_token():
    pass


if __name__ == '__main__':
    """
    serve a flask app returning the closest Hogs and their
    prots. to be integrated into another webapp.

    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-names", help="nameglob for all profilers to be served", type =str)

    args = vars(parser.parse_args(sys.argv[1:]))

    lsh_forests = glob.glob(  config_utils.datadir + args['names'] )
    forests = { :profiler( ) }

    app.run(debug=True)
