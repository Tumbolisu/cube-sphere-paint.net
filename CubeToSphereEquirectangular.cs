// Name: Cube -> Sphere (Equirectangular)
// Submenu: Projection
// Author: Tumby#5171
// Title: Project Cubemap to Spheremap
// Version: 1.0
// Desc: Turns a cubemap into an equirectangular spheremap, eliminating edge-distortions, but introducing pole-distortions.
// Keywords: projection|cube|cubemap|sphere|spheremap|equirectangular|skybox
// URL:
// Help: This effect uses the six images of a cube to generate a texture that can be directly used on a UV-sphere.        The six images must be placed as follows:        (Ultra short description:) [up,dn,--,--], [rt,ft,lf,bk].        (Long description:) In the top left corner is the "up" image. To the right of it lies the "down" image. The rest of the top row is unused. Below the "up" image is the "right" image. To the right are the "front", "left" and "back" images, in order. The "up" and "right" images should line up nicely. The "right", "front", "left" and "back" images should also line up perfectly.        Note 1: Each image should be a perfect square, which also means that the full layout should be twice as wide as it is tall.        Note 2: If you were to place the "down" image below the "right" image, it would line up and create the net of the cube. However, such a layout takes up too much space.        Note 3: The names of the 6 images are based off of Source engine skyboxes, for which this effect was made.
#region UICode
IntSliderControl user_yaw_offset = 0; // [-180,180] Rotate Yaw (Left/Right)
IntSliderControl user_samples = 5; // [1,32] Super-Sampling Size (1 = Disable)
ListBoxControl user_interpolation_choice = 1; // Interpolation Type|Nearest Neighbour (Crisp, Aliased)|Bilinear (Blurry, Antialiased)
ListBoxControl user_window_choice = 1; // Super-Sampling Window Type|Box (Simple, Blurry)|Sinc (Sharper)
#endregion

/*******************************************************************************
    A Paint.NET plugin which converts cubemaps to equirectangular spheremaps.
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

// The following variables never change, so they are calculated once in the Pre-Render.

// required for bilinear sampling.
Surface up, dn, rt, ft, lf, bk;
//      +Z, -Z, +Y, +X, -Y, -X

// super-sampling window/filter.
double[,] window;

void PreRender(Surface dst, Surface src)
{
    int w = src.Width;  // full image width
    int h = src.Height;  // full image height
    int tile_width = w / 4;  // texture/tile width
    int tile_height = h / 2;  // texture/tile height

    // Tiles need to be square, otherwise the following code will be a mess.
    // So, I'm taking the smaller of tile_width and tile_height as the true size.
    int tile_len = (tile_width < tile_height) ? tile_width : tile_height;


    // Work Surfaces //

    Size wrk_size = new Size(tile_len + 2, tile_len + 2);
    int wrk_end = tile_len + 1;

    if (up == null)  up = new Surface(wrk_size);
    if (dn == null)  dn = new Surface(wrk_size);
    if (rt == null)  rt = new Surface(wrk_size);
    if (ft == null)  ft = new Surface(wrk_size);
    if (lf == null)  lf = new Surface(wrk_size);
    if (bk == null)  bk = new Surface(wrk_size);

    // Main Area
    for (int y = 0; y < tile_len; y++)
    {
        for (int x = 0; x < tile_len; x++)
        {
            up[x+1, y+1] = src[0 * tile_len + x, 0 * tile_len + y];
            dn[x+1, y+1] = src[1 * tile_len + x, 0 * tile_len + y];
            rt[x+1, y+1] = src[0 * tile_len + x, 1 * tile_len + y];
            ft[x+1, y+1] = src[1 * tile_len + x, 1 * tile_len + y];
            lf[x+1, y+1] = src[2 * tile_len + x, 1 * tile_len + y];
            bk[x+1, y+1] = src[3 * tile_len + x, 1 * tile_len + y];
        }
    }

    // Edges
    for (int i = 1; i <= tile_len; i++)
    {
        // Left Edge
        up[0, i] = bk[       i,        1];  // Left up = Top bk
        dn[0, i] = bk[       i, tile_len];  // Left dn = Bottom bk
        rt[0, i] = bk[tile_len,        i];  // Left rt = Right bk
        ft[0, i] = rt[tile_len,        i];  // Left ft = Right rt
        lf[0, i] = ft[tile_len,        i];  // Left lf = Right ft
        bk[0, i] = lf[tile_len,        i];  // Left bk = Right lf

        // Right Edge
        up[wrk_end, i] = ft[i,        1];  // Right up = Top ft
        dn[wrk_end, i] = ft[i, tile_len];  // Right dn = Bottom ft
        rt[wrk_end, i] = ft[1,        i];  // Right rt = Left ft
        ft[wrk_end, i] = lf[1,        i];  // Right ft = Left lf
        lf[wrk_end, i] = bk[1,        i];  // Right lf = Left bk
        bk[wrk_end, i] = rt[1,        i];  // Right bk = Left rt

        // Top Edge
        up[i, 0] = lf[       i,        1];  // Top up = Top lf
        dn[i, 0] = rt[       i, tile_len];  // Top dn = Bottom rt
        rt[i, 0] = up[       i, tile_len];  // Top rt = Bottom up
        ft[i, 0] = up[tile_len,        i];  // Top ft = Right up
        lf[i, 0] = up[       i,        1];  // Top lf = Top up
        bk[i, 0] = up[       1,        i];  // Top bk = Left up

        // Bottom Edge
        up[i, wrk_end] = rt[       i,        1];  // Bottom up = Top rt
        dn[i, wrk_end] = lf[       i, tile_len];  // Bottom dn = Bottom lf
        rt[i, wrk_end] = dn[       i,        1];  // Bottom rt = Top dn
        ft[i, wrk_end] = dn[tile_len,        i];  // Bottom ft = Right dn
        lf[i, wrk_end] = dn[       i, tile_len];  // Bottom lf = Bottom dn
        bk[i, wrk_end] = dn[       1,        i];  // Bottom bk = Left dn
    }

    // Corners
    // Not important to "get right", so lets just blend the 2 closest pixels.
    byte alpha = (byte)128;
    int s0,t0,t1,s2;
    for (int i = 0; i < 4; i++)
    {
        switch(i)
        {
        case 0:  // Top-Left Corner
            s0 = 0;
            t0 = 0;
            t1 = t0+1;
            s2 = s0+1;
            break;
        case 1:  // Top-Right Corner
            s0 = wrk_end;
            t0 = 0;
            t1 = t0+1;
            s2 = s0-1;
            break;
        case 2:  // Bottom-Left Corner
            s0 = 0;
            t0 = wrk_end;
            t1 = t0-1;
            s2 = s0+1;
            break;
        default:  // Bottom-Right Corner
            s0 = wrk_end;
            t0 = wrk_end;
            t1 = t0-1;
            s2 = s0-1;
            break;
        }
        
        up[s0,t0] = ColorBgra.Blend(up[s0,t1], up[s2,t0], alpha);
        dn[s0,t0] = ColorBgra.Blend(dn[s0,t1], dn[s2,t0], alpha);
        rt[s0,t0] = ColorBgra.Blend(rt[s0,t1], rt[s2,t0], alpha);
        ft[s0,t0] = ColorBgra.Blend(ft[s0,t1], ft[s2,t0], alpha);
        lf[s0,t0] = ColorBgra.Blend(lf[s0,t1], lf[s2,t0], alpha);
        bk[s0,t0] = ColorBgra.Blend(bk[s0,t1], bk[s2,t0], alpha);
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
        window = null;
        return;
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
        if (up != null)
        {
            up.Dispose();
            up = null;
        }
        
        if (dn != null)
        {
            dn.Dispose();
            dn = null;
        }
        
        if (rt != null)
        {
            rt.Dispose();
            rt = null;
        }
        
        if (ft != null)
        {
            ft.Dispose();
            ft = null;
        }
        
        if (lf != null)
        {
            lf.Dispose();
            lf = null;
        }
        
        if (bk != null)
        {
            bk.Dispose();
            bk = null;
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
    int tile_width = w / 4;  // texture/tile width
    int tile_height = h / 2;  // texture/tile height
    // Tiles need to be square, otherwise code will be a mess.
    // So, I'm taking the smaller of tile_width and tile_height as the true size.
    int tile_len = (tile_width < tile_height) ? tile_width : tile_height;

    double[] vec = new double[3];  // 3D World Vector
    double[] uv = new double[2];  // 2D Image Vector

    double[] pixel_sum;  // Super-sampling sum of pixels. Order: BGRA
    ColorBgra pix = ColorBgra.Black;  // Work Pixel

    char direction;
    
    double y = 0.0;
    double x = 0.0;

    double pitch = 0.0;  // angle made from image y position. ranges from -pi/2 to +pi/2.
    double yaw = 0.0;  // angle made from image x position. ranges from 0 to 2pi. rotations are OK.

    double sin_pitch;
    double cos_pitch;

    ref Surface wrk = ref up;  // This will only ever be a REFERENCE to one of the up,dn,rt,ft,lf,bk surfaces.

    // Some nuisance that has to be explained about the conversion from pixel-coordinate to angle:
    // The X-coordinates range from 0 to w-1.
    // If you image each pixel to be a square, (for the purpose of super-sampling, they literally are)
    // then the ACTUAL coordinates extend in every direction by up to 0.5 pixels.
    // So, instead of mapping the values from (0, w-1) to (0, 2*PI), we use (-0.5, w-0.5) to (0, 2*PI).
    // The same thing happens with the Y-coordinates.


    for (int yy = rect.Top; yy < rect.Bottom; yy++)
    {
        if (IsCancelRequested) return;

        for (int xx = rect.Left; xx < rect.Right; xx++)
        {
            pixel_sum = new double[4];

            for (int sample_y = 0; sample_y < user_samples; sample_y++)
            {
                y = yy + Remap(sample_y, -0.5, user_samples-0.5, -0.5, 0.5);
                pitch = Remap(y, -0.5, h-0.5, 0.5*Math.PI, -0.5*Math.PI);
                sin_pitch = Math.Sin(pitch);
                cos_pitch = Math.Cos(pitch);

                for (int sample_x = 0; sample_x < user_samples; sample_x++)
                {
                    x = xx + Remap(sample_x, -0.5, user_samples-0.5, -0.5, 0.5);
                    yaw = Remap(x, -0.5, w-0.5, 0, 2*Math.PI);
                    yaw -= (user_yaw_offset + 45.0) * (Math.PI / 180.0);

                    vec[0] = cos_pitch * Math.Sin(yaw);
                    vec[1] = cos_pitch * Math.Cos(yaw);
                    vec[2] = sin_pitch;


                    if (Math.Abs(vec[0]) > Math.Abs(vec[1]))  // |X| > |Y|  (X or Z is biggest)
                    {
                        if (Math.Abs(vec[0]) > Math.Abs(vec[2]))  // |X| > |Z|  (X is biggest)
                        {
                            direction = (vec[0] >= 0.0) ? 'X' : 'x';
                        }
                        else  // (Z is biggest)
                        {
                            direction = (vec[2] >= 0.0) ? 'Z' : 'z';
                        }
                    }
                    else  // |Y| >= |X|  (Y or Z biggest)
                    {
                        if (Math.Abs(vec[1]) > Math.Abs(vec[2]))  // |Y| > |Z|  (Y is biggest)
                        {
                            direction = (vec[1] >= 0.0) ? 'Y' : 'y';
                        }
                        else  // (Z is biggest)
                        {
                            direction = (vec[2] >= 0.0) ? 'Z' : 'z';
                        }
                    }
					

					switch (direction)
					{
                    case 'X':
						uv[0] = - vec[1] / vec[0];
						uv[1] = - vec[2] / vec[0];
                        wrk = ref ft;
						break;
						
                    case 'x':
						uv[0] = - vec[1] / vec[0];
						uv[1] =   vec[2] / vec[0];
                        wrk = ref bk;
						break;
						
                    case 'Y':
						uv[0] =   vec[0] / vec[1];
						uv[1] = - vec[2] / vec[1];
                        wrk = ref rt;
						break;
						
                    case 'y':
						uv[0] =   vec[0] / vec[1];
						uv[1] =   vec[2] / vec[1];
                        wrk = ref lf;
						break;
						
                    case 'Z':
						uv[0] =   vec[0] / vec[2];
						uv[1] =   vec[1] / vec[2];
                        wrk = ref up;
						break;
						
                    case 'z':
						uv[0] = - vec[0] / vec[2];
						uv[1] =   vec[1] / vec[2];
                        wrk = ref dn;
						break;

                    default:
                        return;
					}

                    
                    // Remap from World-Vector to Texture-Vector.
                    // uv[0] = Remap(uv[0], -1.0, 1.0, 0, tile_len-1);
                    // uv[1] = Remap(uv[1], -1.0, 1.0, 0, tile_len-1);
                    uv[0] = Remap(uv[0], -1.0, 1.0, -0.5, tile_len-0.5);
                    uv[1] = Remap(uv[1], -1.0, 1.0, -0.5, tile_len-0.5);

                    // offset for surface sampling
                    uv[0] += 1.0;
                    uv[1] += 1.0;

                    switch(user_interpolation_choice)
                    {
                    case 0:  // Nearest-Neighbour
                        pix = wrk[(int)Math.Round(uv[0]), (int)Math.Round(uv[1])];
                        break;

                    case 1:  // Bilinear
                        pix = wrk.GetBilinearSampleClamped((float)uv[0], (float)uv[1]);
                        break;

                    default:  // Invalid Interpolation!
                        return;
                    }

                    pixel_sum[0] += pix.B * window[sample_x, sample_y];
                    pixel_sum[1] += pix.G * window[sample_x, sample_y];
                    pixel_sum[2] += pix.R * window[sample_x, sample_y];
                    pixel_sum[3] += pix.A * window[sample_x, sample_y];
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
