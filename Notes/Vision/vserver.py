#!/usr/local/bin/python 

import numpy as np
import cv2
from flask import Flask, render_template, Response

app = Flask(__name__)

def captureFrame():
    print("Capturing 1")
    cap = cv2.VideoCapture(-1)
    if ( cap.isOpened()):
        print("Cameral is open()")
    else:
        print("Camera Failed to open")
        return b'0'

    cap.set(3, 640)
    cap.set(4, 480)
    trials = 10
    frame = np.array([0])
    stat  = False
    print("Capturing 3")
    
    #while ( not stat and trials > 0 and (frame == 0).all() ):
    while ( trials > 0 ):
        stat, frame = cap.read()
        trials -= 1;
    
    cap.release()

    ret = b'0'
    if (stat):
        ret = cv2.imencode('.jpg', frame)[1].tobytes()
    
    return ret;

def captureCam(width=640, height=480, iterations=None):
    cap = cv2.VideoCapture(0)
    cap.set(3,width) # set Width
    cap.set(4,height) # set Height
    while(True):
        ret, frame = cap.read()
        #frame = cv2.flip(frame, -1) # Flip camera vertically
        #gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

        if ( (frame == 0).all() ):
            print("Zeros ...")
        cv2.imshow('frame', frame)
        #cv2.imshow('gray', gray)

        k = cv2.waitKey(30) & 0xff
        if k == 27: # press 'ESC' to quit
            break
    cap.release()
    cv2.destroyAllWindows()

#captureCam(640,480)
@app.route('/')
def index():
    jpg = captureFrame()
    return Response(jpg, mimetype='image/jpeg')

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=False, threaded=False)
