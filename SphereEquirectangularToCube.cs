// Name: Sphere (Equirectangular) -> Cube
// Submenu: Projection
// Author: Tumby#5171
// Title: Project Spheremap to Cubemap
// Version: 2.0
// Desc: Turns an equirectangular spheremap into a cubemap, eliminating pole-distortions, but introducing edge-distortions.
// Keywords: projection|cube|cubemap|sphere|spheremap|equirectangular|skybox
// URL: [See RTF File]
// Help: [See RTF File]
#region UICode
IntSliderControl user_yaw_offset = 0; // [-180,180] Rotate Yaw (Left/Right)
IntSliderControl user_samples = 5; // [1,32] Super-Sampling Size (1 = Disable)
ListBoxControl user_window_choice = 1; // Super-Sampling Window Type|Box (Simple, Blurry)|Sinc (Sharper)
ListBoxControl user_interpolation_choice = 1; // Interpolation Type|Nearest Neighbour (Crisp, Aliased)|Bilinear (Blurry, Antialiased)
CheckboxControl user_hdr = false; // Compressed HDR
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

// Blending compressed HDR colors is quite different. This function takes care of that.
ColorBgra ColorBlendHDR(ColorBgra pix0, ColorBgra pix1, double frac)
{
    if (frac <= 0) { return pix0; }
    if (frac >= 1) { return pix1; }

    double B = (pix0.B * pix0.A) * (1 - frac) + (pix1.B * pix1.A) * frac;
    double G = (pix0.G * pix0.A) * (1 - frac) + (pix1.G * pix1.A) * frac;
    double R = (pix0.R * pix0.A) * (1 - frac) + (pix1.R * pix1.A) * frac;

    double M = (B > G ? (B > R ? B : R) : (G > R ? G : R));

    pix0.A = (byte)Math.Ceiling(M / 255);

    pix0.R = (byte)(R / pix0.A);
    pix0.G = (byte)(G / pix0.A);
    pix0.B = (byte)(B / pix0.A);

    return pix0;
}

ColorBgra BilinearSampleClamped(Surface surf, double x, double y)
{
    if (user_hdr)
    {
        x = Math.Clamp(x, 0, surf.Width-1);
        y = Math.Clamp(y, 0, surf.Height-1);

        int florx = (int)Math.Floor(x);
        int ceilx = (int)Math.Ceiling(x);
        int flory = (int)Math.Floor(y);
        int ceily = (int)Math.Ceiling(y);

        ColorBgra pix0 = surf[florx, flory];
        ColorBgra pix1 = surf[florx, ceily];
        ColorBgra pix2 = surf[ceilx, flory];
        ColorBgra pix3 = surf[ceilx, ceily];
        
        x = x - florx;
        y = y - flory;

        pix0 = ColorBlendHDR(pix0, pix1, y);
        pix2 = ColorBlendHDR(pix2, pix3, y);

        return ColorBlendHDR(pix0, pix2, x);
    }
    else
    {
        return surf.GetBilinearSampleClamped((float)x, (float)y);
    }
}

// For Debugging.
ColorBgra ColorHDRToLDR(ColorBgra pix, byte scale)
{
    int t;

    t = (pix.R * pix.A) / scale;
    pix.R = (byte)(t > 255 ? 255 : t);

    t = (pix.G * pix.A) / scale;
    pix.G = (byte)(t > 255 ? 255 : t);

    t = (pix.B * pix.A) / scale;
    pix.B = (byte)(t > 255 ? 255 : t);

    pix.A = (byte)255;

    return pix;
}

// Helper Class for super-sampling pixels
class SuperSample
{
    private double B,G,R,A;
    private bool user_hdr;

    public SuperSample(bool is_hdr)
    {
        B = G = R = A = 0;
        user_hdr = is_hdr;
    }

    public void AddPixel(ColorBgra pix, double weight)
    {
        if (user_hdr)
        {
            B += pix.B * pix.A * weight;
            G += pix.G * pix.A * weight;
            R += pix.R * pix.A * weight;
        }
        else
        {
            B += pix.B * weight;
            G += pix.G * weight;
            R += pix.R * weight;
            A += pix.A * weight;
        }
    }

    public ColorBgra ToColorBgra()
    {
        ColorBgra pix = new ColorBgra();
        if (user_hdr)
        {
            double M = (B > G ? (B > R ? B : R) : (G > R ? G : R));
            pix.A = (byte)Math.Ceiling(M / 255);
            pix.R = (byte)(R / pix.A);
            pix.G = (byte)(G / pix.A);
            pix.B = (byte)(B / pix.A);
        }
        else
        {
            pix.B = (byte)Math.Round(B);
            pix.G = (byte)Math.Round(G);
            pix.R = (byte)Math.Round(R);
            pix.A = (byte)Math.Round(A);
        }
        return pix;
    }
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

    SuperSample super;  // Super-sampling sum of pixels.
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
            super = new SuperSample(user_hdr);

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
                        pix = BilinearSampleClamped(wrk, s, t);
                        break;

                    default:  // Invalid Interpolation!
                        return;
                    }

                    super.AddPixel(pix, window[sample_x, sample_y]);
                }
            }  // end of super-sampling

            pix = super.ToColorBgra();
            dst[xx,yy] = pix;
        }
    }  // end of pixel loops
}