// Name: Sphere (Equirectangular) -> Cube
// Submenu: Projection
// Author: Tumby#5171
// Title: Project Spheremap to Cubemap
// Version: 1.0
// Desc: Turns an equirectangular spheremap into a cubemap, eliminating pole-distortions, but introducing edge-distortions.
// Keywords: projection|cube|cubemap|sphere|spheremap|equirectangular|skybox
// URL:
// Help: This effect takes a spheremap, which is a texture that can be directly applied to a UV-sphere, and splits it into the six images of a corresponding cubemap.        The six images will be placed as follows:        (Ultra short description:) [up,dn,--,--], [rt,ft,lf,bk].        (Long description:) In the top left corner is the "up" image. To the right of it lies the "down" image. The rest of the top row is unused. Below the "up" image is the "right" image. To the right are the "front", "left" and "back" images, in order. The "up" and "right" images should line up nicely. The "right", "front", "left" and "back" images should also line up perfectly.        Note 1: In order for each image to be a square, the input image must be twice as wide as it is tall. Using other ratios might cause problems, for which I have no responsibility.        Note 2: If you were to place the "down" image below the "right" image, it would line up and create the net of the cube. However, such a layout takes up too much space.        Note 3: The names of the 6 images are based off of Source engine skyboxes, for which this effect was made.
#region UICode
IntSliderControl user_yaw_offset = 0; // [-180,180] Rotate Yaw (Left/Right)
IntSliderControl user_samples = 5; // [1,32] Super-Sampling Size (1 = Disable)
ListBoxControl user_interpolation_choice = 1; // Interpolation Type|Nearest Neighbour (Crisp, Aliased)|Bilinear (Blurry, Antialiased)
ListBoxControl user_window_choice = 1; // Super-Sampling Window Type|Box (Simple, Blurry)|Sinc (Sharper)
#endregion

/*******************************************************************************
    A Paint.NET plugin which converts equirectangular spheremaps to cubemaps.
    Copyright (C) 2022  R.B. aka "Tumby" aka "Tumbolisu"

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
	
	You can contact me via tumbolisu@gmx.de.
*******************************************************************************/

// Utility function.
// Map a value from one range to another linearly.
double Remap(double value, double in_min, double in_max, double out_min, double out_max)
{
    return out_min + (((value - in_min) * (out_max - out_min)) / (in_max - in_min));
}

// Utility function.
// I literally have to implement fmod from scratch, unless I make this a Visual Studio based plugin.
double FMod(double dividend, double divisor)
{
    return dividend - Math.Truncate(dividend / divisor) * divisor;
}

double FModPositive(double dividend, double divisor)
{
    double result = dividend - Math.Truncate(dividend / divisor) * divisor;
    if (result < 0.0) result += divisor;
    return result;
}


// The following variables never change, so they are calculated once in the Pre-Render.

// A work surface is required for bilinear sampling.
// It's simply a copy of the src surface with a 1-pixel border added.
// It simply repeats in the x direction.
// The y direction needs special care.
Surface wrk;

// super-sampling window/filter.
double[,] window;

void PreRender(Surface dst, Surface src)
{
    // Work Surface //

    Size wrk_size = src.Size;
    wrk_size.Height += 2;
    wrk_size.Width += 2;

    if (wrk == null)
    {
        wrk = new Surface(wrk_size);
    }

    // Main Area
    for (int y = 0; y < src.Height; y++)
    {
        for (int x = 0; x < src.Width; x++)
        {
            wrk[x+1,y+1] = src[x,y];
        }
    }

    // Left and Right Columns
    for (int y = 0; y < src.Height; y++)
    {
        wrk[0, y+1] = src[src.Width-1, y];
        wrk[wrk.Width-1, y+1] = src[0, y];
    }

    // Top and Bottom Rows
    for (int x = 0; x < wrk.Width; x++)
    {
        int x_rhs = (x + (src.Width / 2)) % src.Width;  // 180 degree shift
        wrk[x, 0] = wrk[x_rhs, 1];
        wrk[x, wrk.Height-1] = wrk[x_rhs, wrk.Height-2];
    }


    // Sampling Window //

    // Make a new Array no matter what. Otherwise an old array might stick around after changing user_samples, which causes OutOfBounds.
    window = new double[user_samples, user_samples];

    double[] window_1d = new double[user_samples];

    switch(user_window_choice)
    {
    case 0:  // Box
        for (int i = 0; i < user_samples; i++)
        {
            window_1d[i] = 1.0 / user_samples;
        }
        break;

    case 1:  // Sinc
        {
            double t;
            double sum = 0.0;
            double eps = 1.0 / 4096.0;
            for (int i = 0; i < user_samples; i++)
            {
                t = Remap(i, -0.5, user_samples-0.5, -Math.PI, Math.PI);
                if (t > eps || t < -eps)
                {
                    window_1d[i] = Math.Sin(t) / t;
                }
                else
                {
                    window_1d[i] = 1.0;
                }
                sum += window_1d[i];
            }
            for (int i = 0; i < user_samples; i++)
            {
                window_1d[i] /= sum;
            }
        }
        break;

    default:  // Invalid Window Type!
        {
            window = null;
            return;
        }
    }

    for (int i = 0; i < user_samples; i++)
    {
        for (int j = 0; j < user_samples; j++)
        {
            window[i,j] = window_1d[i] * window_1d[j];
        }
    }
}

protected override void OnDispose(bool disposing)
{
    if (disposing)
    {
        // Release any surfaces or effects you've created.
        if (wrk != null)
        {
            wrk.Dispose();
            wrk = null;
        }
    }
    
    base.OnDispose(disposing);
}

void Render(Surface dst, Surface src, Rectangle rect)
{
    // These 4 variables never change, and could be done globally with the Pre-Render.
    // But just think of how much time would be wasted from the shared thread access.
    // They are literally just 4 ints.
    int w = src.Width;  // full image width
    int h = src.Height;  // full image height
    int tile_w = w / 4;  // texture/tile width
    int tile_h = h / 2;  // texture/tile height

    int[] tile_offset = new int[2];  // which texture to sample

    double[] vec = new double[3];  // 3D World Vector
    double[] uv = new double[2];  // 2D Image Vector

    double[] pixel_sum;  // Super-sampling sum of pixels. Order: BGRA
    ColorBgra pix = ColorBgra.Black;  // Work Pixel

    double pitch = 0.0;  // angle made from image y position. ranges from -pi/2 to +pi/2.
    double yaw = 0.0;  // angle made from image x position. ranges from 0 to 2pi. rotations are OK.

    double y = 0.0;
    double x = 0.0;

    for (int yy = rect.Top; yy < rect.Bottom; yy++)
    {
        if (IsCancelRequested) return;

        for (int xx = rect.Left; xx < rect.Right; xx++)
        {
            pixel_sum = new double[4];

            for (int sample_y = 0; sample_y < user_samples; sample_y++)
            {
                y = yy + Remap(sample_y, -0.5, user_samples-0.5, -0.5, 0.5);

                for (int sample_x = 0; sample_x < user_samples; sample_x++)
                {
                    x = xx + Remap(sample_x, -0.5, user_samples-0.5, -0.5, 0.5);
                    
                    // Find which tile this pixel is on.
                    // This is necessary to get the components of vec right.

                    tile_offset[0] = 0;
                    tile_offset[1] = 0;

                    for (int i = 1; i < 4; i++)
                    {
                        if (xx >= tile_w * i)
                        {
                            tile_offset[0] += 1;
                        }
                    }

                    for (int i = 1; i < 2; i++)
                    {
                        if (yy >= tile_h * i)
                        {
                            tile_offset[1] += 1;
                        }
                    }

                    // Texture-Space Coordinate.
                    uv[0] = x - tile_offset[0] * tile_w;
                    uv[1] = y - tile_offset[1] * tile_h;
                    
                    // Preparation for building the World-Space Vector.
                    // uv[0] = Remap(uv[0], 0, tile_w-1, -1.0, 1.0);
                    // uv[1] = Remap(uv[1], 0, tile_h-1, -1.0, 1.0);
                    uv[0] = Remap(uv[0], -0.5, tile_w-0.5, -1.0, 1.0);
                    uv[1] = Remap(uv[1], -0.5, tile_h-0.5, -1.0, 1.0);

                    // Create a point on the unit-cube.

                    if (tile_offset[1] == 0)  // top row
                    {
                        switch(tile_offset[0])
                        {
                        case 0:  // +Z, up
                            vec[0] =  uv[0];
                            vec[1] =  uv[1];
                            vec[2] =  1.0;
                            break;

                        case 1:  // -Z, down
                            vec[0] =  uv[0];
                            vec[1] = -uv[1];
                            vec[2] = -1.0;
                            break;

                        default:  // pixel is outside the texture layout
                            pix = ColorBgra.Transparent;
                            dst[xx,yy] = pix;
                            continue;
                        }
                    }
                    else  // bottom row
                    {
                        switch(tile_offset[0])
                        {
                        case 0:  // +Y
                            vec[0] =  uv[0];
                            vec[1] =  1.0;
                            vec[2] = -uv[1];
                            break;

                        case 1:  // +X
                            vec[0] =  1.0;
                            vec[1] = -uv[0];
                            vec[2] = -uv[1];
                            break;

                        case 2:  // -Y
                            vec[0] = -uv[0];
                            vec[1] = -1.0;
                            vec[2] = -uv[1];
                            break;

                        case 3:  // -X
                            vec[0] = -1.0;
                            vec[1] =  uv[0];
                            vec[2] = -uv[1];
                            break;
                        }
                    }

                    // Normalize vec to create a point on the unit-sphere.
                    double length = Math.Sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
                    vec[0] /= length;
                    vec[1] /= length;
                    vec[2] /= length;

                    // Now we can calculate the pitch and yaw.
                    // vec[0] = Cos(pitch) * Sin(yaw);
                    // vec[1] = Cos(pitch) * Cos(yaw);
                    // vec[2] = Sin(pitch);
                    
                    // yaw is calculated from the X and Y components by using atan2.
                    yaw = Math.Atan2(vec[0], vec[1]);  // -pi to +pi  (360 deg)
                    yaw += (user_yaw_offset + 45) * (Math.PI / 180.0);
                    yaw = FModPositive(yaw, 2*Math.PI);

                    // pitch is calculated pretty easily from the Z component.
                    pitch = Math.Asin(vec[2]);  // -pi/2 to +pi/2  (180 deg)

                    double s = Remap(yaw, 0, 2*Math.PI, -0.5, w-0.5);
                    //s = FModPositive(s, w);  This isn't harmful, but it's unnecessary.

                    double t = Remap(pitch, 0.5*Math.PI, -0.5*Math.PI, -0.5, h-0.5);
                    //t = FModPositive(t, h);  This causes extreme angles to wrap around the picture, which is NOT WANTED for the pitch!

                    // offset for sampling from work surface
                    s += 1.0;
                    t += 1.0;

                    switch(user_interpolation_choice)
                    {
                    case 0:  // Nearest-Neighbour
                        pix = wrk[(int)Math.Round(s), (int)Math.Round(t)];
                        break;

                    case 1:  // Bilinear
                        pix = wrk.GetBilinearSampleClamped((float)s, (float)t);
                        break;

                    default:  // Invalid Interpolation!
                        return;
                    }

                    pixel_sum[0] += pix.B * window[sample_x, sample_y];
                    pixel_sum[1] += pix.G * window[sample_x, sample_y];
                    pixel_sum[2] += pix.R * window[sample_x, sample_y];
                    pixel_sum[3] += pix.A * window[sample_x, sample_y];

                    //pix.R = (byte)(s*16);
                    //pix.G = (byte)0;
                    //pix.G = (byte)(t*32);
                    //pix.B = (s < 0 || s > 255 || t < 0 || t > 255) ? (byte)255 : (byte)0;
                }
            }  // end of super-sampling

            
            pix.B = (byte)Math.Round(pixel_sum[0]);
            pix.G = (byte)Math.Round(pixel_sum[1]);
            pix.R = (byte)Math.Round(pixel_sum[2]);
            pix.A = (byte)Math.Round(pixel_sum[3]);
            

            dst[xx,yy] = pix;
        }
    }  // end of pixel loops
}