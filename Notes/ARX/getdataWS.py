#!/anaconda/bin/python
from flask import Flask, jsonify
import json

from flask_cors import CORS
app = Flask(__name__)
CORS(app)

from influxdb import InfluxDBClient
client = InfluxDBClient(host='localhost',port=8086,username=None,password=None,ssl=False,database='mydb')

@app.route('/', methods=['GET'])
def index():
    return "Version 1.0"


@app.route('/q/<q>/')
def getQ(q):
    print(f'Query {q}')
    
    res=client.query(q)
    ret = res.raw['series'] #[0]['values'][0:1000]
    ret1= json.dumps(ret)
    return ret1;


@app.route('/data/<table>/')
@app.route('/data/<table>/<start>')
@app.route('/data/<table>/<start>/<end>')
def getdata(table, start=None, end=None):
    print(f'Getting Data from {table} {start}-{end}')
    
    res=client.query(f"select a from {table}")
    ret = res.raw['series'][0]['values'][0:1000]
    ret1= json.dumps(ret)
    return ret1;
    

if __name__ == '__main__':
    app.run(debug=True)