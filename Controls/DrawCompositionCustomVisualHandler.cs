using System;
using System.Numerics;
using Avalonia;
using Avalonia.Controls.Shapes;
using Avalonia.Media;
using Avalonia.Platform;
using Avalonia.Rendering.Composition;
using Avalonia.Skia;
using SkiaSharp;

namespace EffectsDemo;

internal class DrawCompositionCustomVisualHandler : CompositionCustomVisualHandler
{
    private bool _running;
    private Stretch? _stretch;
    private StretchDirection? _stretchDirection;
    private Size? _size;
    private readonly object _sync = new();
    private SKPaint? _paint;
    private SKShader? _shader;
    private SKRuntimeEffectUniforms? _uniforms;
    private SKRuntimeEffect? _effect;
    private SKRuntimeEffectChildren? _children;
    private float _time;

    public bool CanDraw = false;

    public override void OnMessage(object message)
    {
        if (message is not DrawPayload msg)
        {
            return;
        }

        switch (msg)
        {
            case
            {
                HandlerCommand: HandlerCommand.Start,
                Animation: { } an,
                Size: { } size,
                Stretch: { } st,
                StretchDirection: { } sd
            }:
            {
                _running = true;
                _size = size;
                _stretch = st;
                _stretchDirection = sd;
                RegisterForNextAnimationFrameUpdate();
                break;
            }
            case
            {
                HandlerCommand: HandlerCommand.Update,
                Size: { } size,
                Stretch: { } st,
                StretchDirection: { } sd
            }:
            {
                _size = size;
                _stretch = st;
                _stretchDirection = sd;
                RegisterForNextAnimationFrameUpdate();
                break;
            }
            case
            {
                HandlerCommand: HandlerCommand.Stop
            }:
            {
                _running = false;
                break;
            }
            case
            {
                HandlerCommand: HandlerCommand.Dispose
            }:
            {
                DisposeImpl();
                break;
            }
        }
    }

    public override void OnAnimationFrameUpdate()
    {
        if (!_running)
            return;

        Invalidate();
        RegisterForNextAnimationFrameUpdate();
    }

    private void DisposeImpl()
    {
        lock (_sync)
        {
            // TODO:
        }
    }

    private void CreatePaint()
    {
        using var image = SKImage.FromEncodedData(AssetLoader.Open(new Uri("avares://EffectsDemoDRM/Assets/mandrill.png")));
        using var imageShader = image.ToShader();

// https://www.shadertoy.com/view/tdlXDM#
        var shaderCode5 =
"""
uniform vec3      iResolution;           // viewport resolution (in pixels)
uniform float     iTime;                 // shader playback time (in seconds)
uniform float     iTimeDelta;            // render time (in seconds)
uniform float     iFrameRate;            // shader frame rate
uniform int       iFrame;                // shader playback frame
uniform float     iChannelTime[4];       // channel playback time (in seconds)
uniform vec3      iChannelResolution[4]; // channel resolution (in pixels)
uniform vec4      iMouse;                // mouse pixel coords. xy: current (if MLB down), zw: click
//uniform samplerXX iChannel0..3;          // input channel. XX = 2D/Cube
uniform vec4      iDate;                 // (year, month, day, time in seconds)
uniform float     iSampleRate;           // sound sample rate (i.e., 44100)
// uniform float     PI;

////////////////////////////////////////////////////////////////////////////////
//
// Retro cube - mimic a vectors-cube on an old CRT-display
//
// Copyright 2019 Mirco Müller
//
// Author(s):
//   Mirco "MacSlow" Müller <macslow@gmail.com>
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3, as published
// by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranties of
// MERCHANTABILITY, SATISFACTORY QUALITY, or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

const vec3 red = vec3 (1.0, 0.0, 0.0);
const vec3 green = vec3 (0.0, 1.0, 0.0);
const vec3 blue = vec3 (0.0, 0.0, 1.0);
const vec3 white = vec3 (1.0);
const vec3 orange = vec3 (1.0, 0.4, 0.125);
const vec3 black = vec3 (0.2, 0.3, 0.2);
const vec3 cyan = vec3 (0.0, 1.0, 1.0);
const vec3 magenta = vec3 (1.0, 0.0, 1.0);
const vec3 yellow = vec3 (1.0, 1.0, 0.0);
const float SIZE = .003;

float distLine (vec2 p, vec2 a, vec2 b) {
    vec2 pa = p - a;
    vec2 ba = b - a;
    float t = clamp ( dot (pa, ba) / dot (ba, ba), .0, 1.);
    return length (pa - ba*t);
}

float lineMask (vec2 uv, vec2 a, vec2 b) {
    float d = distLine (uv, a, b);
    float thickness = SIZE;
    return smoothstep (thickness, .125*thickness, d);
}

vec3 glowLine (vec2 uv, vec2 a, vec2 b, vec3 rgbGlow) {
    float m = lineMask (uv, a, b);
    float dist = distLine (uv, a, b);
    float brightness = SIZE/pow (.085 + 2.*dist, 2.);
    vec3 color = m*vec3 (.7);
    color += rgbGlow*brightness;
	return color;
}

struct boxType {vec4 p[8];};

mat4 trans (vec3 t)
{
    mat4 mat = mat4 (vec4 (1., .0, .0, .0),
                     vec4 (.0, 1., .0, .0),
                     vec4 (.0, .0, 1., .0),
                     vec4 (t.x, t.y, t.z, 1.));
    return mat;
}

mat4 rotX (float angle)
{
    float rad = radians (angle);
    float c = cos (rad);
    float s = sin (rad);

    mat4 mat = mat4 (vec4 (1., .0, .0, .0),
                     vec4 (.0,   c,   s, .0),
                     vec4 (.0,  -s,   c, .0),
                     vec4 (.0, .0, .0, 1.));

    return mat;
}

mat4 rotY (float angle)
{
    float rad = radians (angle);
    float c = cos (rad);
    float s = sin (rad);

    mat4 mat = mat4 (vec4 (  c, .0,  -s, .0),
                     vec4 (.0, 1., .0, .0),
                     vec4 (  s, .0,   c, .0),
                     vec4 (.0, .0, .0, 1.));

    return mat;
}

mat4 rotZ (float angle)
{
    float rad = radians (angle);
    float c = cos (rad);
    float s = sin (rad);

    mat4 mat = mat4 (vec4 (  c,   s, .0, .0),
                     vec4 ( -s,   c, .0, .0),
                     vec4 (.0, .0, 1.0, .0),
                     vec4 (.0, .0, .0, 1.));

    return mat;
}

vec4 main (vec2 fragCoord) {
  	vec4 fragColor;
    vec2 uv = (fragCoord/iResolution.xy)* 2. - 1.;
    uv.x *= iResolution.x/iResolution.y;
    uv *= 1. + .4*length(uv);
    uv.x += .0035*cos(40.*uv.y + 2.*iTime);

    boxType box;
    box.p[0] = vec4 ( 0.1,  0.1,  0.1, 1.0);
    box.p[1] = vec4 ( 0.1, -0.1,  0.1, 1.0);
    box.p[2] = vec4 (-0.1, -0.1,  0.1, 1.0);
    box.p[3] = vec4 (-0.1,  0.1,  0.1, 1.0);
    box.p[4] = vec4 ( 0.1,  0.1, -0.1, 1.0);
    box.p[5] = vec4 ( 0.1, -0.1, -0.1, 1.0);
    box.p[6] = vec4 (-0.1, -0.1, -0.1, 1.0);
    box.p[7] = vec4 (-0.1,  0.1, -0.1, 1.0);

    float t = 8. + 14.*iTime;
    mat4 rot3d = rotX (-4.*t)*rotY (3.*t)*rotZ (2.*t);
    mat4 model = trans (vec3 (.0, .0, -.275))*rot3d;
    box.p[0] = model * box.p[0];
    box.p[1] = model * box.p[1];
    box.p[2] = model * box.p[2];
    box.p[3] = model * box.p[3];
    box.p[4] = model * box.p[4];
    box.p[5] = model * box.p[5];
    box.p[6] = model * box.p[6];
    box.p[7] = model * box.p[7];

    vec3 boxCol = glowLine (uv, box.p[0].xy / box.p[0].z, box.p[1].xy / box.p[1].z, red);
    boxCol += glowLine (uv, box.p[1].xy / box.p[1].z, box.p[2].xy / box.p[2].z, green);
    boxCol += glowLine (uv, box.p[2].xy / box.p[2].z, box.p[3].xy / box.p[3].z, orange);
    boxCol += glowLine (uv, box.p[3].xy / box.p[3].z, box.p[0].xy / box.p[0].z, cyan);
    boxCol += glowLine (uv, box.p[4].xy / box.p[4].z, box.p[5].xy / box.p[5].z, blue);
    boxCol += glowLine (uv, box.p[5].xy / box.p[5].z, box.p[6].xy / box.p[6].z, red);
    boxCol += glowLine (uv, box.p[6].xy / box.p[6].z, box.p[7].xy / box.p[7].z, yellow);
    boxCol += glowLine (uv, box.p[7].xy / box.p[7].z, box.p[4].xy / box.p[4].z, green);
    boxCol += glowLine (uv, box.p[0].xy / box.p[0].z, box.p[4].xy / box.p[4].z, blue);
    boxCol += glowLine (uv, box.p[1].xy / box.p[1].z, box.p[5].xy / box.p[5].z, cyan);
    boxCol += glowLine (uv, box.p[2].xy / box.p[2].z, box.p[6].xy / box.p[6].z, green);
    boxCol += glowLine (uv, box.p[3].xy / box.p[3].z, box.p[7].xy / box.p[7].z, magenta);

    boxCol = boxCol / (1. + boxCol);
    boxCol = sqrt (boxCol);
    boxCol *= mix (1., .5, .5 + .5*cos (500.*uv.y));

    fragColor = vec4(boxCol, 1.);

    return fragColor;
}
""";

var shaderCode =
"""
uniform vec3      iResolution;           // viewport resolution (in pixels)
uniform float     iTime;                 // shader playback time (in seconds)
uniform float     iTimeDelta;            // render time (in seconds)
uniform float     iFrameRate;            // shader frame rate
uniform int       iFrame;                // shader playback frame
uniform float     iChannelTime[4];       // channel playback time (in seconds)
uniform vec3      iChannelResolution[4]; // channel resolution (in pixels)
uniform vec4      iMouse;                // mouse pixel coords. xy: current (if MLB down), zw: click
//uniform samplerXX iChannel0..3;          // input channel. XX = 2D/Cube
uniform vec4      iDate;                 // (year, month, day, time in seconds)
uniform float     iSampleRate;           // sound sample rate (i.e., 44100)
// uniform float     PI;

const int SLICES = 4;
const float DENSITY = 3.0;

vec3 saturation(vec3 c, float t) {
    return mix(vec3(dot(c,vec3(0.2126,0.7152,0.0722))),c,t);
}
mat3 fromEuler(vec3 ang) {
	vec2 a1 = vec2(sin(ang.x),cos(ang.x));
    vec2 a2 = vec2(sin(ang.y),cos(ang.y));
    vec2 a3 = vec2(sin(ang.z),cos(ang.z));
    mat3 m;
    m[0] = vec3(a1.y*a3.y+a1.x*a2.x*a3.x,a1.y*a2.x*a3.x+a3.y*a1.x,-a2.y*a3.x);
	m[1] = vec3(-a2.y*a1.x,a1.y*a2.y,a2.x);
	m[2] = vec3(a3.y*a1.x*a2.x+a1.y*a3.x,a1.x*a3.x-a1.y*a3.y*a2.x,a2.y*a3.y);
	return m;
}
bool intersectionRayBox(vec3 o, vec3 d, vec3 ext, out vec3 r0, out vec3 r1) {
    vec3 t0 = (-o - ext) / d;
    vec3 t1 = (-o + ext) / d;
    vec3 n = min(t0,t1); n.x = max(max(n.x,n.y),n.z);
    vec3 f = max(t0,t1); f.x = min(min(f.x,f.y),f.z);
    r0 = o + d * n.x;
    r1 = o + d * f.x;
    return bool(step(n.x,f.x));
}

float integrationFunc(float x, float a) {
    return x * 0.5 - cos(x * a) / (2.0 * a);
}

float functionMean(float a, float b, float f) {
    a = -a * 0.5 + 0.5;
    b = -b * 0.5 + 0.5;
    float Fa = integrationFunc(a,f);
    float Fb = integrationFunc(b,f);
    return (Fb - Fa) / (b - a);
}

// main
vec4 main( vec2 fragCoord ) {
    vec2 uv = (-iResolution.xy + 2.0 * fragCoord) / iResolution.y;
    vec2 mouse = iMouse.xy / iResolution.xy * 4.0 - 2.0;

    // ray
    vec3 ang;
    if(iMouse.z > 0.0) {
        ang = vec3(0.0,-mouse.y,mouse.x);
    } else {
        ang = vec3(sin(iTime*0.4)*2.0,0.0,cos(iTime*0.35)*3.0);
    }
	mat3 rot = fromEuler(ang);
    vec3 ori = vec3(0.0,0.0,4.5) * rot;
    vec3 dir = normalize(vec3(uv.xy,-3.0)) * rot;

    // color
    vec3 p, rp0, rp1;
    vec3 rcolor = vec3(0.0);
    if(intersectionRayBox(ori,dir,vec3(1.0),rp0,rp1)) {
        for(int i = 0; i < SLICES; i++) {
            vec3 r0 = rp0;
            vec3 r1 = rp1;

            r0 = mix(r0,r1,float(SLICES-i-1)/float(SLICES));
            r1 = r0 + (rp1-rp0)/float(SLICES);
            r0 += 1.7;
            r1 += 1.7;

            // modulate color
            float fm = 0.6;
            vec3 color;
            color.x = functionMean(r0.x,r1.x,7.0*fm);
            color.y = functionMean(r0.x,r1.x,11.0*fm);
            color.z = functionMean(r0.x,r1.x,13.0*fm);
            color.yz *= functionMean(r0.x,r1.x,27.0*fm);
            color = 1.0 - color;

            color.z *= functionMean(r0.y,r1.y,11.0*fm);
            color.y *= functionMean(r0.y,r1.y,13.0*fm);
            color.x *= functionMean(r0.y,r1.y,17.0*fm);
            color = 1.0 - color;

            color.z *= functionMean(r0.z,r1.z,5.0*fm);
            color.y *= functionMean(r0.z,r1.z,7.0*fm);
            color.x *= functionMean(r0.z,r1.z,11.0*fm);
            color = 1.0 - color;
            color = pow(color,vec3(8.0));
            color += 0.02;
            color = log(1.0+color*1.5);

            // modulate density
            float d = length(r1 - r0) * functionMean(r0.z+r0.y,r1.z+r1.y,22.0);
            rcolor = mix(rcolor, color, clamp(log(1.0+d*DENSITY),0.0,1.0));
        }

        // sRGB
        rcolor = pow(rcolor, vec3(1.0/2.2));
        rcolor = saturation(rcolor, 2.0);
    }

	return vec4(rcolor,1.0);

}

""";

var shaderCode3 =
"""
uniform vec3      iResolution;           // viewport resolution (in pixels)
uniform float     iTime;                 // shader playback time (in seconds)
uniform float     iTimeDelta;            // render time (in seconds)
uniform float     iFrameRate;            // shader frame rate
uniform int       iFrame;                // shader playback frame
uniform float     iChannelTime[4];       // channel playback time (in seconds)
uniform vec3      iChannelResolution[4]; // channel resolution (in pixels)
uniform vec4      iMouse;                // mouse pixel coords. xy: current (if MLB down), zw: click
//uniform samplerXX iChannel0..3;          // input channel. XX = 2D/Cube
uniform vec4      iDate;                 // (year, month, day, time in seconds)
uniform float     iSampleRate;           // sound sample rate (i.e., 44100)
// uniform float     PI;

// Copyright (c) 2013 Andrew Baldwin (baldand)
// License = Attribution-NonCommercial-ShareAlike (http://creativecommons.org/licenses/by-nc-sa/3.0/deed.en_US)

// "Mirror Cube"
// A simple ray tracer and a simple scene - one cube


const vec3 up = vec3(0.,1.,0.);

float intersectfloor(vec3 ro, vec3 rd, float height, out float t0)
{
	if (rd.y==0.0) {
		t0 = 100000.0;
		return 0.0;
	}

	t0 = -(ro.y - height)/rd.y;
	t0 = min(100000.0,t0);
	return t0;
}

float intersectbox(vec3 ro, vec3 rd, float size, out float t0, out float t1, out vec3 normal)
// Calculate intersections with origin-centred axis-aligned cube with sides length size
// Returns positive value if there are intersections
{
    vec3 ir = 1.0/rd;
    vec3 tb = ir * (vec3(-size*.5)-ro);
    vec3 tt = ir * (vec3(size*.5)-ro);
    vec3 tn = min(tt, tb);
    vec3 tx = max(tt, tb);
    vec2 t = max(tn.xx, tn.yz);
    t0 = max(t.x, t.y);
    t = min(tx.xx, tx.yz);
    t1 = min(t.x, t.y);
	float d = (t1-t0);
	vec3 i = ro + t0*rd;
	normal = step(size*.499,abs(i))*sign(i);
	if (t0<-0.01) d = t0;
	return d;
}

float intersect(vec3 boxPos, vec3 ro, vec3 rd, out vec3 intersection, out vec3 normal, out int material, out float t)
{
	float tb0=0.0;
	float tb1=0.0;
	vec3 boxnormal;
	float dbox = intersectbox(ro-boxPos,rd,1.,tb0,tb1,boxnormal);
	float tf = 0.0;
	float dfloor = intersectfloor(ro,rd,0.,tf);
	t = tf;
	float d = dfloor;
	material = 0; // Sky
	if (d>=0.) {
		normal = vec3(0.,1.,0.);
		material = 2; // Floor
	}
	if (dbox>=0.) {
		t = tb0;
		d = dbox;
		normal = boxnormal;
		material = 1; // Box
		if (t<0.) d=-0.1;
	}
	intersection = ro+t*rd;
	return d;
}

vec4 main( vec2 fragCoord )
{
	float rotspeed = iTime*1.+iMouse.x/iResolution.x;
	vec3 light = vec3(5.,4.+3.*sin(-rotspeed*.4),2.);
	float radius = sin(rotspeed*.1)*2.+4.;
	vec3 boxPos = vec3(0.3,1.5+sin(rotspeed),0.2);
	vec3 eye = vec3(radius*sin(rotspeed),2.*sin(.1*rotspeed)+2.5+2.*iMouse.y/iResolution.y,radius*cos(rotspeed*1.));
	vec3 screen = vec3((radius-1.)*sin(rotspeed),1.5*sin(.1*rotspeed)+2.+2.*iMouse.y/iResolution.y,(radius-1.)*cos(rotspeed*1.));
    vec2 screenSize = vec2(iResolution.x/iResolution.y,1.0);
	vec2 uv = fragCoord.xy / iResolution.xy;
	vec2 offset = screenSize * (uv - 0.5);
	vec3 right = cross(up,normalize(screen - eye));
	vec3 ro = screen + offset.y*up + offset.x*right;
	vec3 rd = normalize(ro - eye);
	vec3 i = vec3(0.);
	vec3 n = vec3(0.);
	int m,m2;
	float d,lightd,ra,global,direct,shade,t,tlight;
	vec3 lrd,i2,n2;
	vec3 c=vec3(0.);
	vec3 ca=vec3(0.);
	float lra=1.;
	for (int reflections=0;reflections<10;reflections++) {
		// Find the direct ray hit
		d = intersect(boxPos,ro,rd,i,n,m,t);
		// Check for shadows to the light
		lrd = normalize(light-i);
		tlight = length(light-i);
		lightd = smoothstep(.5*length(i-i2),.0,intersect(boxPos,i,lrd,i2,n2,m2,t));
		if (t>tlight) lightd=1.0;
		// Colouring
		global = .3;
		direct = max( (10./length(lrd)) * dot( lrd, n) ,0.0);
		shade = global + direct*lightd;
		if (m==0) { ra=0.0; c = vec3(0.9,2.0,2.5); }
		if (m==1) { ra=0.2; c = shade*(.5+.5*(i-boxPos)); }
		if (m==2) {
			ra = 0.3;
			vec2 mxz = abs(fract(i.xz)*2.-1.);
			float fade = clamp(1.-length(i.xz)*.05,0.,1.);
			float fc =mix(.5,smoothstep(1.,.9,mxz.x+mxz.y),fade);
			c = vec3(fc*shade);
		}
		// Calculate any reflection on the next iteration
		ca += lra*c;
		lra *= ra;
		rd = reflect(rd,n);
		ro = i+0.01*rd;
	}
	return vec4(ca/(1.+ca),1.);
}
""";

var shaderCode4 =
"""
uniform vec3      iResolution;           // viewport resolution (in pixels)
uniform float     iTime;                 // shader playback time (in seconds)
uniform float     iTimeDelta;            // render time (in seconds)
uniform float     iFrameRate;            // shader frame rate
uniform int       iFrame;                // shader playback frame
uniform float     iChannelTime[4];       // channel playback time (in seconds)
uniform vec3      iChannelResolution[4]; // channel resolution (in pixels)
uniform vec4      iMouse;                // mouse pixel coords. xy: current (if MLB down), zw: click
//uniform samplerXX iChannel0..3;          // input channel. XX = 2D/Cube
uniform vec4      iDate;                 // (year, month, day, time in seconds)
uniform float     iSampleRate;           // sound sample rate (i.e., 44100)
// uniform float     PI;

// uncomment for a cross section view
//#define CROSS_SECTION

//------------------------------------------------------------------------
// Camera
//
// Move the camera. In this case it's using time and the mouse position
// to orbitate the camera around the origin of the world (0,0,0), where
// the yellow sphere is.
//------------------------------------------------------------------------
void doCamera( out vec3 camPos, out vec3 camTar, in float time, in float mouseX )
{

    float an = 0.3*iTime + 10.0*mouseX;

	camPos = vec3(4.5*sin(an),2.0,4.5*cos(an));
    camTar = vec3(0.0,0.0,0.0);
}


//------------------------------------------------------------------------
// Background
//
// The background color. In this case it's just a black color.
//------------------------------------------------------------------------
vec3 doBackground( void )
{
    return vec3( 0.0, 0.0, 0.0);
}

// all three basic bodies are symmetric across the XYZ planes
// octahedron and rhombic dodecahedron have been scaled to align
// with the vertices of the cube.

// 1D distance of X Y Z planes
vec2 cube(vec3 p, float r) {
    vec3 o = abs(p);
	float s = o.x;
	s = max(s, o.y);
	s = max(s, o.z);
	return vec2(s-r, 0.0);
}

// 3D distance of XYZ cross diagonal plane
vec2 octahedron(vec3 p, float r) {
    vec3 o = abs(p) / sqrt(3.0);
	float s = o.x+o.y+o.z;
	return vec2(s-r*2.0/sqrt(3.0), 1.0);
}

// 2D distance of XY YZ ZX diagonal planes
vec2 rhombic(vec3 p, float r) {
    vec3 o = abs(p) / sqrt(2.0);
	float s = o.x+o.y;
	s = max(s, o.y+o.z);
	s = max(s, o.z+o.x);
	return vec2(
        s-r*sqrt(2.0),
        2.0);
}

vec2 min2(vec2 a, vec2 b) {
    return (a.x <= b.x)?a:b;
}

vec2 max2(vec2 a, vec2 b) {
    return (a.x > b.x)?a:b;
}

const float SHAPE_COUNT = 8.0;
vec3 get_factors(int i) {
    if (i == 0) {
        // cube
        return vec3(1.0, 6.0/4.0, 1.0);
    } else if (i == 1) {
        // truncated cube
        return vec3(1.0, 6.0/5.0, 1.0);
    } else if (i == 2) {
        // cuboctahedron
        return vec3(1.0, 1.0, 1.0);
    } else if (i == 3) {
        // truncated octahedron
        return vec3(4.0/3.0, 1.0, 1.0);
    } else if (i == 4) {
        // truncated cuboctahedron
        return vec3(sqrt(3.0/2.0), 2.0/sqrt(3.0), 1.0);
    } else if (i == 5) {
        // rhombicuboctahedron
        return vec3(sqrt(2.0), sqrt(5.0/3.0), 1.0);
    } else if (i == 6) {
        // octahedron
        return vec3(2.0, 1.0, 1.0);
    }
    return vec3(0.0);
}

vec2 plane( vec3 p) {
    return vec2(p.y+2.0,3.0);
}

//------------------------------------------------------------------------
// Modelling
//
// Defines the shapes (a sphere in this case) through a distance field, in
// this case it's a sphere of radius 1.
//------------------------------------------------------------------------
vec2 add_plane(vec3 p, vec2 m) {
    return min2(plane(p),m);
}

vec2 doModel( vec3 p ) {
    float k = iTime*0.5;
    //k = 1.0;
    float u = smoothstep(0.0,1.0,smoothstep(0.0,1.0,fract(k)));
    int s1 = int(mod(k,SHAPE_COUNT));
    int s2 = int(mod(k+1.0,SHAPE_COUNT));
    if (s1 == 6) {
        return add_plane(p, mix(octahedron(p, 1.0), rhombic(p, 1.0), u));
    } else if (s1 == 7) {
        return add_plane(p, mix(rhombic(p, 1.0), cube(p, 1.0), u));
    } else {
        vec3 f = mix(get_factors(s1),
                   get_factors(s2), u);
        return add_plane(p, max2(max2(cube(p,f.x),octahedron(p, f.y)), rhombic(p, f.z)));
    }

}

//------------------------------------------------------------------------
// Material
//
// Defines the material (colors, shading, pattern, texturing) of the model
// at every point based on its position and normal. In this case, it simply
// returns a constant yellow color.
//------------------------------------------------------------------------
vec3 doMaterial( in vec3 pos, in vec3 nor )
{
    float k = doModel(pos).y;
    return mix(mix(mix(vec3(1.0,0.07,0.01),vec3(0.2,1.0,0.01),clamp(k,0.0,1.0)),
               vec3(0.1,0.07,1.0),
               clamp(k-1.0,0.0,1.0)),
               vec3(0.1),
               clamp(k-2.0,0.0,1.0));
}

//------------------------------------------------------------------------
// Lighting
//------------------------------------------------------------------------
float calcSoftshadow( in vec3 ro, in vec3 rd );

vec3 doLighting( in vec3 pos, in vec3 nor, in vec3 rd, in float dis, in vec3 mal )
{
    vec3 lin = vec3(0.0);

    // key light
    //-----------------------------
    vec3  lig = normalize(vec3(1.0,0.7,0.9));
    float dif = dot(nor,lig) * 0.5 + 0.5;
    float sha = 0.0; if( dif>0.01 ) sha=calcSoftshadow( pos+0.01*nor, lig );
    lin += dif;


    // surface-light interacion
    //-----------------------------
    vec3 col = mal*lin;


    // fog
    //-----------------------------
	col *= exp(-0.01*dis*dis);

    return col;
}

float calcIntersection( in vec3 ro, in vec3 rd )
{
	const float maxd = 20.0;           // max trace distance
	const float precis = 0.001;        // precission of the intersection
    float h = precis*2.0;
    float t = 0.0;
	float res = -1.0;
    for( int i=0; i<90; i++ )          // max number of raymarching iterations is 90
    {
        if( h<precis||t>maxd ) break;
	    h = doModel( ro+rd*t ).x;
        t += h;
    }

    if( t<maxd ) res = t;
    return res;
}

vec3 calcNormal( in vec3 pos )
{
    const float eps = 0.002;             // precision of the normal computation

    const vec3 v1 = vec3( 1.0,-1.0,-1.0);
    const vec3 v2 = vec3(-1.0,-1.0, 1.0);
    const vec3 v3 = vec3(-1.0, 1.0,-1.0);
    const vec3 v4 = vec3( 1.0, 1.0, 1.0);

	return normalize( v1*doModel( pos + v1*eps ).x +
					  v2*doModel( pos + v2*eps ).x +
					  v3*doModel( pos + v3*eps ).x +
					  v4*doModel( pos + v4*eps ).x );
}

float calcSoftshadow( in vec3 ro, in vec3 rd )
{
    float res = 1.0;
    float t = 0.0005;                 // selfintersection avoidance distance
	float h = 1.0;
    for( int i=0; i<40; i++ )         // 40 is the max numnber of raymarching steps
    {
        h = doModel(ro + rd*t).x;
        res = min( res, 64.0*h/t );   // 64 is the hardness of the shadows
		t += clamp( h, 0.02, 2.0 );   // limit the max and min stepping distances
    }
    return clamp(res,0.0,1.0);
}

mat3 calcLookAtMatrix( in vec3 ro, in vec3 ta, in float roll )
{
    vec3 ww = normalize( ta - ro );
    vec3 uu = normalize( cross(ww,vec3(sin(roll),cos(roll),0.0) ) );
    vec3 vv = normalize( cross(uu,ww));
    return mat3( uu, vv, ww );
}

vec4 main( vec2 fragCoord )
{
    vec2 p = (-iResolution.xy + 2.0*fragCoord.xy)/iResolution.y;
    vec2 m = iMouse.xy/iResolution.xy;

    //-----------------------------------------------------
    // camera
    //-----------------------------------------------------

    // camera movement
    vec3 ro, ta;
    doCamera( ro, ta, iTime, m.x );

    // camera matrix
    mat3 camMat = calcLookAtMatrix( ro, ta, 0.0 );  // 0.0 is the camera roll

	// create view ray
	vec3 rd = normalize( camMat * vec3(p.xy,2.0) ); // 2.0 is the lens length

    //-----------------------------------------------------
	// render
    //-----------------------------------------------------

	vec3 col = doBackground();

	// raymarch
    float t = calcIntersection( ro, rd );
    if( t>-0.5 )
    {
        // geometry
        vec3 pos = ro + t*rd;
        vec3 nor = calcNormal(pos);

        // materials
        vec3 mal = doMaterial( pos, nor );

        col = doLighting( pos, nor, rd, t, mal );
	}

	//-----------------------------------------------------
	// postprocessing
    //-----------------------------------------------------
    // gamma
	col = pow( clamp(col,0.0,1.0), vec3(0.4545) );

    return vec4( col, 1.0 );
}

""";


var src = """
                  uniform vec3      iResolution;           // viewport resolution (in pixels)
                  uniform float     iTime;                 // shader playback time (in seconds)
                  uniform float     iTimeDelta;            // render time (in seconds)
                  uniform float     iFrameRate;            // shader frame rate
                  uniform int       iFrame;                // shader playback frame
                  uniform float     iChannelTime[4];       // channel playback time (in seconds)
                  uniform vec3      iChannelResolution[4]; // channel resolution (in pixels)
                  uniform vec4      iMouse;                // mouse pixel coords. xy: current (if MLB down), zw: click
                  //uniform samplerXX iChannel0..3;          // input channel. XX = 2D/Cube
                  uniform vec4      iDate;                 // (year, month, day, time in seconds)
                  uniform float     iSampleRate;           // sound sample rate (i.e., 44100)

                  // Created by greenbird10
                  // License Creative Commons Attribution-NonCommercial-ShareAlike 3.0

                  float hash(vec2 p) {
                  	return 0.5*(
                      sin(dot(p, vec2(271.319, 413.975)) + 1217.13*p.x*p.y)
                      ) + 0.5;
                  }

                  float noise(vec2 p) {
                    vec2 w = fract(p);
                    w = w * w * (3.0 - 2.0*w);
                    p = floor(p);
                    return mix(
                      mix(hash(p+vec2(0,0)), hash(p+vec2(1,0)), w.x),
                      mix(hash(p+vec2(0,1)), hash(p+vec2(1,1)), w.x), w.y);
                  }

                  // wave octave inspiration
                  // Alexander Alekseev - Seascape
                  // https://www.shadertoy.com/view/Ms2SD1
                  float map_octave(vec2 uv) {
                    uv = (uv + noise(uv)) / 2.5;
                    uv = vec2(uv.x*0.6-uv.y*0.8, uv.x*0.8+uv.y*0.6);
                    vec2 uvsin = 1.0 - abs(sin(uv));
                    vec2 uvcos = abs(cos(uv));
                    uv = mix(uvsin, uvcos, uvsin);
                    float val = 1.0 - pow(uv.x * uv.y, 0.65);
                    return val;
                  }

                  float map(vec3 p) {
                    vec2 uv = p.xz + iTime/2.;
                    float amp = 0.6, freq = 2.0, val = 0.0;
                    for(int i = 0; i < 3; ++i) {
                      val += map_octave(uv) * amp;
                      amp *= 0.3;
                      uv *= freq;
                      // uv = vec2(uv.x*0.6-uv.y*0.8, uv.x*0.8+uv.y*0.6);
                    }
                    uv = p.xz - 1000. - iTime/2.;
                    amp = 0.6, freq = 2.0;
                    for(int i = 0; i < 3; ++i) {
                      val += map_octave(uv) * amp;
                      amp *= 0.3;
                      uv *= freq;
                      // uv = vec2(uv.x*0.6-uv.y*0.8, uv.x*0.8+uv.y*0.6);
                    }
                    return val + 3.0 - p.y;
                  }

                  vec3 getNormal(vec3 p) {
                    float eps = 1./iResolution.x;
                    vec3 px = p + vec3(eps, 0, 0);
                    vec3 pz = p + vec3(0, 0, eps);
                    return normalize(vec3(map(px),eps,map(pz)));
                  }

                  // raymarch inspiration
                  // Alexander Alekseev - Seascape
                  // https://www.shadertoy.com/view/Ms2SD1
                  float raymarch(vec3 ro, vec3 rd, out vec3 outP, out float outT) {
                      float l = 0., r = 26.;
                      const int steps = 16;
                      float dist = 1000000.;
                      for(int i = 0; i < steps; ++i) {
                          float mid = (r+l)/2.;
                          float mapmid = map(ro + rd*mid);
                          dist = min(dist, abs(mapmid));
                          if(mapmid > 0.) {
                          	l = mid;
                          }
                          else {
                          	r = mid;
                          }
                          if(r - l < 1./iResolution.x) break;
                      }
                      outP = ro + rd*l;
                      outT = l;
                      return dist;
                  }

                  float fbm(vec2 n) {
                  	float total = 0.0, amplitude = 1.0;
                  	for (int i = 0; i < 5; i++) {
                  		total += noise(n) * amplitude;
                  		n += n;
                  		amplitude *= 0.4;
                  	}
                  	return total;
                  }

                  float lightShafts(vec2 st) {
                      float angle = -0.2;
                      vec2 _st = st;
                      float t = iTime / 16.;
                      st = vec2(st.x * cos(angle) - st.y * sin(angle),
                                st.x * sin(angle) + st.y * cos(angle));
                      float val = fbm(vec2(st.x*2. + 200. + t, st.y/4.));
                      val += fbm(vec2(st.x*2. + 200. - t, st.y/4.));
                      val = val / 3.;
                      float mask = pow(clamp(1.0 - abs(_st.y-0.15), 0., 1.)*0.49 + 0.5, 2.0);
                      mask *= clamp(1.0 - abs(_st.x+0.2), 0., 1.) * 0.49 + 0.5;
                  	return pow(val*mask, 2.0);
                  }

                  vec2 bubble(vec2 uv, float scale) {
                      if(uv.y > 0.2) return vec2(0.);
                      float t = iTime/4.;
                      vec2 st = uv * scale;
                      vec2 _st = floor(st);
                      vec2 bias = vec2(0., 4. * sin(_st.x*128. + t));
                      float mask = smoothstep(0.1, 0.2, -cos(_st.x*128. + t));
                      st += bias;
                      vec2 _st_ = floor(st);
                      st = fract(st);
                      float size = noise(_st_)*0.07+0.01;
                      vec2 pos = vec2(noise(vec2(t, _st_.y*64.1)) * 0.8 + 0.1, 0.5);
                      if(length(st.xy - pos) < size) {
                          return (st + pos) * vec2(.1, .2) * mask;
                      }
                      return vec2(0.);
                  }

                  //void mainImage(out vec4 fragColor, in vec2 fragCoord) {
                  half4 main(vec2 fragCoord) {
                      vec3 ro = vec3(0.,0.,2.);
                      vec3 lightPos = vec3(8, 3, -3);
                      vec3 lightDir = normalize(lightPos - ro);

                      // adjust uv
                      vec2 uv = fragCoord;
                      uv = (-iResolution.xy + 2.0*uv) / iResolution.y;
                      uv.y *= 0.5;
                      uv.x *= 0.45;
                      uv += bubble(uv, 12.) + bubble(uv, 24.); // add bubbles

                      vec3 rd = normalize(vec3(uv, -1.));
                      vec3 hitPos;
                      float hitT;
                      vec3 seaColor = vec3(11,82,142)/255.;
                      vec3 color;

                      // waves
                      float dist = raymarch(ro, rd, hitPos, hitT);
                      float diffuse = dot(getNormal(hitPos), rd) * 0.5 + 0.5;
                      color = mix(seaColor, vec3(15,120,152)/255., diffuse);
                      color += pow(diffuse, 12.0);
                  	// refraction
                      vec3 ref = normalize(refract(hitPos-lightPos, getNormal(hitPos), 0.05));
                      float refraction = clamp(dot(ref, rd), 0., 1.0);
                      color += vec3(245,250,220)/255. * 0.6 * pow(refraction, 1.5);

                      vec3 col = vec3(0.);
                      col = mix(color, seaColor, pow(clamp(0., 1., dist), 0.2)); // glow edge
                      col += vec3(225,230,200)/255. * lightShafts(uv); // light shafts

                      // tone map
                      col = (col*col + sin(col))/vec3(1.8, 1.8, 1.9);

                      // vignette
                      // inigo quilez - Stop Motion Fox
                      // https://www.shadertoy.com/view/3dXGWB
                      vec2 q = fragCoord / iResolution.xy;
                      col *= 0.7+0.3*pow(16.0*q.x*q.y*(1.0-q.x)*(1.0-q.y),0.2);

                      //fragColor = vec4(col,1.0);
                      return vec4(col,1.0);
                  }
                  """;

        _time = 0f;

        _effect = SKRuntimeEffect.CreateShader(shaderCode, out var errorText);
        _uniforms = new SKRuntimeEffectUniforms(_effect)
        {
            ["iResolution"] = new [] { 512f, 512f, 0f },
            ["iTime"] = _time,
            ["iMouse"] = new [] { 0f, 0f, -1f, -1f },
        };


        //_children = new SKRuntimeEffectChildren(_effect) { ["iImage1"] = imageShader };
        _children = new SKRuntimeEffectChildren(_effect);

        _shader = _effect.ToShader(_uniforms, _children);
        _paint = new SKPaint { Shader = _shader };
    }

    private void UpdatePaint()
    {
        if (_uniforms is { } && _effect is { } && _paint is { })
        {
            _time += 1f / 60f;
            _uniforms["iTime"] = _time;
            _shader?.Dispose();
            _shader = _effect.ToShader(_uniforms, _children);
            _paint.Shader = _shader;
        }
    }

    private void Draw(SKCanvas canvas)
    {
        if (CanDraw) {
            canvas.Save();

            if (_paint is null)
            {
                CreatePaint();
            }
            else
            {
                UpdatePaint();
            }

            // TODO:
            canvas.DrawRect(0, 0, 512, 512, new SKPaint { Color = SKColors.Black });
            canvas.DrawRect(0, 0, 512, 512, _paint);

            canvas.Restore();
        }
    }

    public override void OnRender(ImmediateDrawingContext context)
    {
        lock (_sync)
        {

            if (_stretch is not { } st
                || _stretchDirection is not { } sd)
            {
                return;
            }

            var leaseFeature = context.TryGetFeature<ISkiaSharpApiLeaseFeature>();
            if (leaseFeature is null)
            {
                return;
            }

            var rb = GetRenderBounds();

            var size = _size ?? rb.Size;

            var viewPort = new Rect(rb.Size);
            var sourceSize = new Size(size.Width, size.Height);
            if (sourceSize.Width <= 0 || sourceSize.Height <= 0)
            {
                return;
            }

            var scale = st.CalculateScaling(rb.Size, sourceSize, sd);
            var scaledSize = sourceSize * scale;
            var destRect = viewPort
                .CenterRect(new Rect(scaledSize))
                .Intersect(viewPort);
            var sourceRect = new Rect(sourceSize)
                .CenterRect(new Rect(destRect.Size / scale));

            var bounds = SKRect.Create(new SKPoint(), new SKSize((float)size.Width, (float)size.Height));
            var scaleMatrix = Matrix.CreateScale(
                destRect.Width / sourceRect.Width,
                destRect.Height / sourceRect.Height);
            var translateMatrix = Matrix.CreateTranslation(
                -sourceRect.X + destRect.X - bounds.Top,
                -sourceRect.Y + destRect.Y - bounds.Left);

            using (context.PushClip(destRect))
            using (context.PushPostTransform(translateMatrix * scaleMatrix))
            {
                using var lease = leaseFeature.Lease();
                var canvas = lease?.SkCanvas;
                if (canvas is null)
                {
                    return;
                }
                Draw(canvas);
            }
        }
    }
}
