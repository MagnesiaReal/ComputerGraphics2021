# ComputerGraphics2021

## Practice 1 diagnosis exam
Just complete a table

## Practice 2 XY proyection
This program output is a XY proyection from the OBJ file which faces are triangles.
### image
![alt text](https://github.com/MagnesiaReal/ComputerGraphics2021/blob/main/tests/2021-12-03-18:17:27.png)

## Practice 3 OBJ_TO_PPM
Perform a single Rotation, Translation and Scale of an object and transform to PPM file.
![alt text](https://github.com/MagnesiaReal/ComputerGraphics2021/blob/main/tests/2021-12-03-18:32:24.png)

## Practice 4 Animation 1
Perform a single rotation, translation and scale, then the program perform autimatic Ry rotation and translation.

Without perspective projection: https://www.youtube.com/watch?v=GiTW3XFR-Sg .

The program update with perspective projection: https://youtu.be/zwxQ31q1uYM .

![alt text](https://github.com/MagnesiaReal/ComputerGraphics2021/blob/main/tests/2021-12-03-18:17:17.png)

## Practice 5 Animation 2
Same as practice 4 but now read a VLF file instead OBJ file, also implmented FULLFILL px array.

## Practice 6 Animation 3
Features:

- Face hidding

- Z buffer

- Light atenuation factor

Video: https://youtu.be/MgJ2_KpBR_Y .

![alt text](https://github.com/MagnesiaReal/ComputerGraphics2021/blob/main/tests/2021-12-03-18:16:57.png)

## Curves_activity1
Its a simple curve in 2 Dimentions following the next ecuations:

x(t) = a<sub>x</sub>t<sup>3</sup> + b<sub>x</sub>t<sup>2</sup> + c<sub>x</sub>t + d<sub>x</sub>

y(t) = a<sub>y</sub>t<sup>3</sup> + b<sub>y</sub>t<sup>2</sup> + c<sub>y</sub>t + d<sub>y</sub>

#### Input
./curves \<ax\> \<bx\> \<cx\> \<dx\> \<ay\> \<by\> \<cy\> \<dy\> \<n-lines\>
  
where \<n-lines\> is the number of lines for draw this curve.

![alt text](https://github.com/MagnesiaReal/ComputerGraphics2021/blob/main/tests/2022-01-04-16:17:17.png)

#### Output

![alt text](https://github.com/MagnesiaReal/ComputerGraphics2021/blob/main/tests/curve2D.png)

## Curves_activity2
Curve in 3 Dimentions with perpective projection, this use the next ecuations:

x(t) = a<sub>x</sub>t<sup>3</sup> + b<sub>x</sub>t<sup>2</sup> + c<sub>x</sub>t + d<sub>x</sub>

y(t) = a<sub>y</sub>t<sup>3</sup> + b<sub>y</sub>t<sup>2</sup> + c<sub>y</sub>t + d<sub>y</sub>

z(t) = a<sub>z</sub>t<sup>3</sup> + b<sub>z</sub>t<sup>2</sup> + c<sub>z</sub>t + d<sub>z</sub>

#### Input

./curves \<ax\> \<bx\> \<cx\> \<dx\> \<ay\> \<by\> \<cy\> \<dy\> \<az\> \<bz\> \<cz\> \<dz\>

all of these are coefficients for x y and z.

![alt text](https://github.com/MagnesiaReal/ComputerGraphics2021/blob/main/tests/2022-01-04-16:33:32.png)

#### Output

![alt text](https://github.com/MagnesiaReal/ComputerGraphics2021/blob/main/tests/curve3D.png)

## Curves_Hermite

The next image explains how it is works:

![alt text](https://github.com/MagnesiaReal/ComputerGraphics2021/blob/main/tests/2022-01-04-16:41:36.png)

#### Input

Hermite curve needs four parameters:

./hermite \<p1x\> \<p1y\> \<p1z\> \<p2x\> \<p2y\> \<p2z\> \<h1x\> \<h1y\> \<h1z\> \<h2x\> \<h2y\> \<h2z\>

this program uses 200 fixed lines for the curve.

![alt text](https://github.com/MagnesiaReal/ComputerGraphics2021/blob/main/tests/2022-01-04-16:26:31.png)

#### Output

![alt text](https://github.com/MagnesiaReal/ComputerGraphics2021/blob/main/tests/hermite_curve.png)

## Curves_Bezier

The next image explains how it is works:

![alt text](https://github.com/MagnesiaReal/ComputerGraphics2021/blob/main/tests/2022-01-04-16:43:26.png)

#### Input

Bezier curve needs four parameters:

./bezier \<p1x\> \<p1y\> \<p1z\> \<p2x\> \<p2y\> \<p2z\> \<b1x\> \<b1y\> \<b1z\> \<b2x\> \<b2y\> \<b2z\>

![alt text](https://github.com/MagnesiaReal/ComputerGraphics2021/blob/main/tests/2022-01-04-16:21:09.png)

#### Output

![alt text](https://github.com/MagnesiaReal/ComputerGraphics2021/blob/main/tests/bezier_curve.png)

## Cuves_Surfaces

The next image is a representacion for each parameter of surface program:

![alt text](https://github.com/MagnesiaReal/ComputerGraphics2021/blob/main/tests/2022-01-04-16:48:07.png)

#### Input

![alt text](https://github.com/MagnesiaReal/ComputerGraphics2021/blob/main/tests/2022-01-04-16:28:33.png)

#### Output

![alt text](https://github.com/MagnesiaReal/ComputerGraphics2021/blob/main/tests/surface.png)



