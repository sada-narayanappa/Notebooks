#!/usr/local/bin/python

'''
This is a Simple Streaming video server from https://blog.miguelgrinberg.com/post/video-streaming-with-flask
you may start the server as: 
[~] gunicorn --worker-class gevent --workers 1 --bind 0.0.0.0:5000 tserver:app

and open a browser with link http://localhost:5000
'''

from flask import Flask, render_template, Response
import time
import cv2
import numpy as np

class Camera(object):
    def __init__(self):
        self.frames = [open(f + '.jpg', 'rb').read() for f in ['1', '2', '3']]

    @staticmethod
    def captureFrame():
        return self.frames[0]
    
    def get_frame(self):
        return self.frames[int(time()) % 3]

class Camera_OCV(object):
    def __init__(self):
        #print("**** Opening Camera ....")
        cap = cv2.VideoCapture(0)
        #print("**** Opening Camera 1....")
        #cap.set(3, 640)
        #cap.set(4, 480)
        self.cap = cap;
        #print("**** Opened Camera ....", self.cap.isOpened(), cap.isOpened())

    @staticmethod
    def captureFrame():
        print("Capturing 1")
        cap = cv2.VideoCapture(-1)
        if ( not cap.isOpened()):
            print("Camera Failed to open")
            return b'0'

        cap.set(3, 640)
        cap.set(4, 480)
        trials = 3
        frame = np.array([0])
        stat  = False

        while ( not stat and trials > 0 and (frame == 0).all() ):
            time.sleep(.3) # Sleep for 200 ms for camera to turn on;
            stat, frame = cap.read()
            trials -= 1;

        cap.release()

        ret = b'0'
        if (stat):
            ret = cv2.imencode('.jpg', frame)[1].tobytes()

        return ret;
    
    def get_frame(self):
        if ( not self.cap.isOpened() ):
            return b'0'
        stat, frame = self.cap.read()
        ret = b'0'
        if (stat):
            ret = cv2.imencode('.jpg', frame)[1].tobytes()
        #print(self.cap, self.cap.isOpened(), stat, ret[0:10])
        return ret;

app = Flask(__name__)

@app.route('/')
def index():
    cam = Camera_OCV()
    jpg = cam.captureFrame()
    return Response(jpg, mimetype='image/jpeg')

def gen(camera):
    while True:
        #time.sleep(.1) # Sleep for 200 ms for camera to turn on;
        frame = camera.get_frame()
        yield (b'--frame\r\n'
               b'Content-Type: image/jpeg\r\n\r\n' + frame + b'\r\n')

@app.route('/video')
@app.route('/video_feed')
def video_feed():
    #cam = Camera()
    cam = Camera_OCV()
    return Response(gen(cam), mimetype='multipart/x-mixed-replace; boundary=frame')


if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=True, threaded=False)