#PPM animation

this program create 240 ppm files that represent 10 seconds of animation with 24fps.
firtly you can perform one rotation, translation and scaling, then the program make animation wich consist
of moving the obj down with Yrotation.

the comand I used to make a video is this 

  ffmpeg -framerate 24 -i "%03d.ppm" -i <songname> -shortest output.avi

-i => inputs 
-shortest => cut the output video when at least one of the inputs finish
-framerate => framerate output video
