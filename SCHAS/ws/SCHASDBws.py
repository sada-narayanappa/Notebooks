#!/usr/local/bin/python

from flask import Flask, jsonify, abort, request, make_response, url_for
from flask_restful import Resource, Api
from json import dumps
import sys, os, importlib, inspect
import re;
import json;
import os
import numpy as np
from IPython.display import HTML
import pandas as pd

PyUtilsPATH="../../PyUtils";
sys.path.append(PyUtilsPATH)
SQL_TXT = "/opt/SCHAS/data/sql/SQL.txt";

from DBUtils import  DBUtils
from GenUtils import *

#db = DBUtils(conn='postgresql://postgres:postgres@localhost:5400/SCHASDB')
db = DBUtils(conn='postgresql://postgres:postgres@localhost:5432/SCHASDB')
db.load(file=SQL_TXT);
#--------------------------------------------------------------------------
js = '''
$("#t1 th").resizable({handles:''})

$("#t1 th").each(function (i) {
    var th = $(this);
    th.css("background-color", '#ededed');
    th.css("text-align", 'center');
});
'''
def _getQueryHTML(q, limit=25, htmlname='dbtable', rType='html' ):
    global js, df;
    
    df = db.execQ(q, limit =limit)
    if ( type(df) != pd.core.frame.DataFrame):
        ret ="Something wrong: executing {} got {}".format(q, df)
        return ret, df
    
    res= df.to_html();
    #res = res.encode('utf-8')
    res=res.replace('\n', '')
    res = res.replace("<table ", "<table id={} width='100%' ".format(htmlname))
    res = res.replace("<td", "<td contenteditable nowrap")
    
    jjs = js.replace("#t1", "#"+htmlname)
    res = res + "<script>" + jjs + "</script>"
    return res, df
    
#--------------------------------------------------------------------------
def getQuery(parm=None):
    req = {'q': "select 'NO QUERY' as HMMM;", 'type': 'html', 'max': '100', 'tname':'dbtable'}
    r = _get(request, parm, req)

    try:
        #print(r)
        ret, df = _getQueryHTML(r.q, int(r.max), r.tname)
        if (r.type.startswith('html') ):
            return str(ret)
        if (True or r.type.startswith('j') ):
            return str(df.to_json());
    except(Exception) as e:
        ret = "Error: " + str(e)
        
    #print(ret)
    return str(ret)

#--------------------------------------------------------------------------
def getQCache(parm=None):
    req = {'id': '1'}
    r = _get(request, parm, req)
    ret = db.QCACHE[r.id]
    return ret
#--------------------------------------------------------------------------
def getTablePK(parm=None):
    req = {'t': 'test'}
    r = _get(request, parm, req)
    try:
        ret = db.insp.get_pk_constraint(r.t)
        print (str(ret))
    except(Exception) as e:
        ret = "Error: " + str(e)
        
    return str(ret)
#--------------------------------------------------------------------------
def getTable(parm=None):
    req = {'t': '', 'type': 'html', 'max': 100, 'tname':'dbtable'}
    r = _get(request, parm, req)
    try:
        #print(r)
        q= "select * from "+r.t 
        ret, df = _getQueryHTML(q, int(r.max), r.tname)
        if (r.type == 'html'):
            return str(ret)
        if (True or r.type.startswith('j') ):
            return str(df.to_json());
    except(Exception) as e:
        ret = "Error: " + str(e)
        
    #print(ret)
    return str(ret)
#--------------------------------------------------------------------------
def getTables():
    e=sorted(db.insp.get_table_names())
    return str(e)
#--------------------------------------------------------------------------
def _get(request, parm, req):
    if (not request):
        return parm or req;
    rreq = parm or req

    for t in request.args.keys():
        rreq[t] = request.args[t]

    return DDict(rreq)
#--------------------------------------------------------------------------
def _default():
    print ("SCHASDB Web service\n\n");
    print (getQuery(1) );

#--------------------------------------------------------------------------
if __name__ == '__main__':
    _default();
