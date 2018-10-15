#!/usr/local/bin/python 
import numpy as np
import cv2
import os

def save(path, image, jpg_quality=95, png_compression=4):
    
    if path.endswith("jpg"):
        cv2.imwrite(path, image, [int(cv2.IMWRITE_JPEG_QUALITY), jpg_quality])
    elif path.endswith(".png"):
        cv2.imwrite(path, image, [int(cv2.IMWRITE_PNG_COMPRESSION), png_compression])
    else:
        cv2.imwrite(path, image)

def captureCam(width=640, height=480, iterations=None):
    cap = cv2.VideoCapture(0)
    cap.set(3,width) # set Width
    cap.set(4,height) # set Height
    while(True and cap.isOpened()):
        ret, frame = cap.read()
        #frame = cv2.flip(frame, -1) # Flip camera vertically
        #gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

        if ( not ret or (frame == 0).all() ):
            print("Zeros ...", ret)
        else:
            file = "test.jpg"
            if ( not os.path.exists(file)):
                save(file, frame)
            
        cv2.imshow('frame', frame)
        #cv2.imshow('gray', gray)

        k = cv2.waitKey(30) & 0xff
        if k == 27: # press 'ESC' to quit
            break
    cap.release()
    cv2.destroyAllWindows()

captureCam(640,480)