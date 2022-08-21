# Cube and Sphere Projections for Paint.NET
A few Paint.NET plugins for converting between images of unwrapped cubes and projected spheres.

<br/>

Includes the following plugins:
### ![CubeToSphereEquirectangular](/CubeToSphereEquirectangular.png "Icon for Cube To Sphere (Equirectangular)") Cube -> Sphere (Equirectangular)
Turns a cubemap into an equirectangular spheremap, eliminating edge-distortions, but introducing pole-distortions.
### ![SphereEquirectangularToCube](/SphereEquirectangularToCube.png "Icon for Sphere (Equirectangular) To Cube") Sphere (Equirectangular) -> Cube
Turns an equirectangular spheremap into a cubemap, eliminating pole-distortions, but introducing edge-distortions.


## Installation
Extract the contents of the ZIP archive to any folder. (It should be a DLL and a BAT file in every ZIP.)

Then simply double-click the BAT file to run it, which will install the DLL file for you.

You can savely remove the ZIP, DLL and BAT files afterwards.


## Usage

All plugins are found in the Effects tab, in the submenu "Projection".

For more information on the usage of a certain plugin, use the plugin and click the ? button in the top-right of the window title bar.


<br/><br/><br/>
I am considering making another set of conversions that use Fisheye Dome projections of a sphere.

This allows for a specific region of the image to have as little distortion as possible. Perfect for painting points of interest, like a sun in the sky.

Of course I won't make this unless there is actually interest in it. (And if I have the time for it as well.)
