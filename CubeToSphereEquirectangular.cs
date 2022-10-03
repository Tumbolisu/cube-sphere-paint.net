// Name: Cube -> Sphere (Equirectangular)
// Submenu: Projection
// Author: Tumby#5171
// Title: Project Cubemap to Spheremap
// Version: 2.0
// Desc: Turns a cubemap into an equirectangular spheremap, eliminating edge-distortions, but introducing pole-distortions.
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

ColorBgra ColorBlendAuto(ColorBgra pix0, ColorBgra pix1, double frac)
{
    if (user_hdr)
    {
        return ColorBlendHDR(pix0, pix1, frac);
    }
    else
    {
        return ColorBgra.Blend(pix0, pix1, (byte)Math.Round(frac*255));
    }
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
    int TL = tile_len;  // Abbreviation. Only use when appropiate!

    // Work Surfaces //

    Size wrk_size = new Size(tile_len + 2, tile_len + 2);
    int wrk_end = tile_len + 1;
    int WE = wrk_end;  // Abbreviation. Only use when appropiate!

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
            up[x+1, y+1] = src[0 * TL + x, 0 * TL + y];
            dn[x+1, y+1] = src[1 * TL + x, 0 * TL + y];
            rt[x+1, y+1] = src[0 * TL + x, 1 * TL + y];
            ft[x+1, y+1] = src[1 * TL + x, 1 * TL + y];
            lf[x+1, y+1] = src[2 * TL + x, 1 * TL + y];
            bk[x+1, y+1] = src[3 * TL + x, 1 * TL + y];
        }
    }

    // Edges
    int mi;  // mirrored direction of i
    for (int i = 1; i <= tile_len; i++)
    {
        mi = tile_len - i;

        // Left Edge
        up[ 0, i] = bk[ i,  1];  // Left up = Top bk
        dn[ 0, i] = bk[mi, TL];  // Left dn = Bottom bk
        rt[ 0, i] = bk[TL,  i];  // Left rt = Right bk
        ft[ 0, i] = rt[TL,  i];  // Left ft = Right rt
        lf[ 0, i] = ft[TL,  i];  // Left lf = Right ft
        bk[ 0, i] = lf[TL,  i];  // Left bk = Right lf

        // Right Edge
        up[WE, i] = ft[mi,  1];  // Right up = Top ft
        dn[WE, i] = ft[ i, TL];  // Right dn = Bottom ft
        rt[WE, i] = ft[ 1,  i];  // Right rt = Left ft
        ft[WE, i] = lf[ 1,  i];  // Right ft = Left lf
        lf[WE, i] = bk[ 1,  i];  // Right lf = Left bk
        bk[WE, i] = rt[ 1,  i];  // Right bk = Left rt

        // Top Edge
        up[ i, 0] = lf[mi,  1];  // Top up = Top lf
        dn[ i, 0] = rt[ i, TL];  // Top dn = Bottom rt
        rt[ i, 0] = up[ i, TL];  // Top rt = Bottom up
        ft[ i, 0] = up[TL, mi];  // Top ft = Right up
        lf[ i, 0] = up[mi,  1];  // Top lf = Top up
        bk[ i, 0] = up[ 1,  i];  // Top bk = Left up

        // Bottom Edge
        up[i, WE] = rt[ i,  1];  // Bottom up = Top rt
        dn[i, WE] = lf[mi, TL];  // Bottom dn = Bottom lf
        rt[i, WE] = dn[ i,  1];  // Bottom rt = Top dn
        ft[i, WE] = dn[TL,  i];  // Bottom ft = Right dn
        lf[i, WE] = dn[mi, TL];  // Bottom lf = Bottom dn
        bk[i, WE] = dn[ 1, mi];  // Bottom bk = Left dn
    }

    // Corners
    // Not important to "get right", so lets just blend the 2 closest pixels.
    double alpha = 0.5;
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
        
        up[s0,t0] = ColorBlendAuto(up[s0,t1], up[s2,t0], alpha);
        dn[s0,t0] = ColorBlendAuto(dn[s0,t1], dn[s2,t0], alpha);
        rt[s0,t0] = ColorBlendAuto(rt[s0,t1], rt[s2,t0], alpha);
        ft[s0,t0] = ColorBlendAuto(ft[s0,t1], ft[s2,t0], alpha);
        lf[s0,t0] = ColorBlendAuto(lf[s0,t1], lf[s2,t0], alpha);
        bk[s0,t0] = ColorBlendAuto(bk[s0,t1], bk[s2,t0], alpha);
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
            double sum = 0;
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
                    window_1d[i] = 1;
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

    SuperSample super;  // Super-sampling sum of pixels.
    ColorBgra pix = ColorBgra.Black;  // Work Pixel

    char direction;
    
    double y = 0;
    double x = 0;

    double pitch = 0;  // angle made from image y position. ranges from -pi/2 to +pi/2.
    double yaw = 0;  // angle made from image x position. ranges from 0 to 2pi. rotations are OK.

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
            super = new SuperSample(user_hdr);

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
                        pix = BilinearSampleClamped(wrk, uv[0], uv[1]);
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
