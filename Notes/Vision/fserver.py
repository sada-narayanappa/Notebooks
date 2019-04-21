#!/usr/local/bin/python

from flask import Flask, render_template, Response
import time
import cv2
import os

def getFile(file, contents):
    ret = None
    if (not os.path.exists(file+".1") and os.path.exists(file)):
        ret = open(file, 'rb').read()
        os.remove(file)
    return ret

app = Flask(__name__)

@app.route('/')
def index():
    return "Video Feed"

def feed(file='test.jpg'):
    contents = None
    while True:
        contents = getFile(file, contents)
        if ( contents ):
            yield (b'--frame\r\n'
                   b'Content-Type: image/jpeg\r\n\r\n' + contents + b'\r\n')
        else:
            time.sleep(3)

@app.route('/video')
def video_feed():
    return Response(feed('test.jpg'), mimetype='multipart/x-mixed-replace; boundary=frame')

if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=True, threaded=False)