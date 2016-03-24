%trainingImageLabeler
trainCascadeObjectDetector('stopSignDetector.xml',positiveInstances,'C:\Users\Dolen\Documents\MATLAB\SrProject\scans\negative','FalseAlarmRate',0.2,'NumCascadeStages',5);
detector = vision.CascadeObjectDetector('stopSignDetector.xml');
img = imread('scans/scan3.jpg');
img = imresize(img,[400 800]);
bbox = step(detector,img);
detectedImg = insertObjectAnnotation(img,'rectangle',bbox,'resistor');
figure;
imshow(detectedImg);
delete('stopSignDetector.xml');