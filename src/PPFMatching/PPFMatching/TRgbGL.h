
#pragma once

#ifndef TRGBGL_H
#define COLORTRGBGL_H

typedef struct TColorGL {
	float red;
	float green;
	float blue;
} TColorGL;

//const int RGBSIZE = 656;
#define RGBSIZE 656

#define COLORaliceBlue           {0.94118f,0.97255f,1.00000f}
#define COLORantiqueWhite        {0.98039f,0.92157f,0.84314f}
#define COLORantiqueWhite1       {1.00000f,0.93725f,0.85882f}
#define COLORantiqueWhite2       {0.93333f,0.87451f,0.80000f}
#define COLORantiqueWhite3       {0.80392f,0.75294f,0.69020f}
#define COLORantiqueWhite4       {0.54510f,0.51373f,0.47059f}
#define COLORaquamarine          {0.49804f,1.00000f,0.83137f}
#define COLORaquamarine1         {0.49804f,1.00000f,0.83137f}
#define COLORaquamarine2         {0.46275f,0.93333f,0.77647f}
#define COLORaquamarine3         {0.40000f,0.80392f,0.66667f}
#define COLORaquamarine4         {0.27059f,0.54510f,0.45490f}
#define COLORazure               {0.94118f,1.00000f,1.00000f}
#define COLORazure1              {0.94118f,1.00000f,1.00000f}
#define COLORazure2              {0.87843f,0.93333f,0.93333f}
#define COLORazure3              {0.75686f,0.80392f,0.80392f}
#define COLORazure4              {0.51373f,0.54510f,0.54510f}
#define COLORbeige               {0.96078f,0.96078f,0.86275f}
#define COLORbisque              {1.00000f,0.89412f,0.76863f}
#define COLORbisque1             {1.00000f,0.89412f,0.76863f}
#define COLORbisque2             {0.93333f,0.83529f,0.71765f}
#define COLORbisque3             {0.80392f,0.71765f,0.61961f}
#define COLORbisque4             {0.54510f,0.49020f,0.41961f}
#define COLORblack               {0.00000f,0.00000f,0.00000f}
#define COLORblanchedAlmond      {1.00000f,0.92157f,0.80392f}
#define COLORblue                {0.00000f,0.00000f,1.00000f}
#define COLORblue1               {0.00000f,0.00000f,1.00000f}
#define COLORblue2               {0.00000f,0.00000f,0.93333f}
#define COLORblue3               {0.00000f,0.00000f,0.80392f}
#define COLORblue4               {0.00000f,0.00000f,0.54510f}
#define COLORblueViolet          {0.54118f,0.16863f,0.88627f}
#define COLORbrown               {0.64706f,0.16471f,0.16471f}
#define COLORbrown1              {1.00000f,0.25098f,0.25098f}
#define COLORbrown2              {0.93333f,0.23137f,0.23137f}
#define COLORbrown3              {0.80392f,0.20000f,0.20000f}
#define COLORbrown4              {0.54510f,0.13725f,0.13725f}
#define COLORburlywood           {0.87059f,0.72157f,0.52941f}
#define COLORburlywood1          {1.00000f,0.82745f,0.60784f}
#define COLORburlywood2          {0.93333f,0.77255f,0.56863f}
#define COLORburlywood3          {0.80392f,0.66667f,0.49020f}
#define COLORburlywood4          {0.54510f,0.45098f,0.33333f}
#define COLORcadetBlue           {0.37255f,0.61961f,0.62745f}
#define COLORcadetBlue1          {0.59608f,0.96078f,1.00000f}
#define COLORcadetBlue2          {0.55686f,0.89804f,0.93333f}
#define COLORcadetBlue3          {0.47843f,0.77255f,0.80392f}
#define COLORcadetBlue4          {0.32549f,0.52549f,0.54510f}
#define COLORchartreuse          {0.49804f,1.00000f,0.00000f}
#define COLORchartreuse1         {0.49804f,1.00000f,0.00000f}
#define COLORchartreuse2         {0.46275f,0.93333f,0.00000f}
#define COLORchartreuse3         {0.40000f,0.80392f,0.00000f}
#define COLORchartreuse4         {0.27059f,0.54510f,0.00000f}
#define COLORchocolate           {0.82353f,0.41176f,0.11765f}
#define COLORchocolate1          {1.00000f,0.49804f,0.14118f}
#define COLORchocolate2          {0.93333f,0.46275f,0.12941f}
#define COLORchocolate3          {0.80392f,0.40000f,0.11373f}
#define COLORchocolate4          {0.54510f,0.27059f,0.07451f}
#define COLORcoral               {1.00000f,0.49804f,0.31373f}
#define COLORcoral1              {1.00000f,0.44706f,0.33725f}
#define COLORcoral2              {0.93333f,0.41569f,0.31373f}
#define COLORcoral3              {0.80392f,0.35686f,0.27059f}
#define COLORcoral4              {0.54510f,0.24314f,0.18431f}
#define COLORcornflowerBlue      {0.39216f,0.58431f,0.92941f}
#define COLORcornsilk            {1.00000f,0.97255f,0.86275f}
#define COLORcornsilk1           {1.00000f,0.97255f,0.86275f}
#define COLORcornsilk2           {0.93333f,0.90980f,0.80392f}
#define COLORcornsilk3           {0.80392f,0.78431f,0.69412f}
#define COLORcornsilk4           {0.54510f,0.53333f,0.47059f}
#define COLORcyan                {0.00000f,1.00000f,1.00000f}
#define COLORcyan1               {0.00000f,1.00000f,1.00000f}
#define COLORcyan2               {0.00000f,0.93333f,0.93333f}
#define COLORcyan3               {0.00000f,0.80392f,0.80392f}
#define COLORcyan4               {0.00000f,0.54510f,0.54510f}
#define COLORdarkBlue            {0.00000f,0.00000f,0.54510f}
#define COLORdarkCyan            {0.00000f,0.54510f,0.54510f}
#define COLORdarkGoldenrod       {0.72157f,0.52549f,0.04314f}
#define COLORdarkGoldenrod1      {1.00000f,0.72549f,0.05882f}
#define COLORdarkGoldenrod2      {0.93333f,0.67843f,0.05490f}
#define COLORdarkGoldenrod3      {0.80392f,0.58431f,0.04706f}
#define COLORdarkGoldenrod4      {0.54510f,0.39608f,0.03137f}
#define COLORdarkGray            {0.66275f,0.66275f,0.66275f}
#define COLORdarkGreen           {0.00000f,0.39216f,0.00000f}
#define COLORdarkGrey            {0.66275f,0.66275f,0.66275f}
#define COLORdarkKhaki           {0.74118f,0.71765f,0.41961f}
#define COLORdarkMagenta         {0.54510f,0.00000f,0.54510f}
#define COLORdarkOliveGreen      {0.33333f,0.41961f,0.18431f}
#define COLORdarkOliveGreen1     {0.79216f,1.00000f,0.43922f}
#define COLORdarkOliveGreen2     {0.73725f,0.93333f,0.40784f}
#define COLORdarkOliveGreen3     {0.63529f,0.80392f,0.35294f}
#define COLORdarkOliveGreen4     {0.43137f,0.54510f,0.23922f}
#define COLORdarkOrange          {1.00000f,0.54902f,0.00000f}
#define COLORdarkOrange1         {1.00000f,0.49804f,0.00000f}
#define COLORdarkOrange2         {0.93333f,0.46275f,0.00000f}
#define COLORdarkOrange3         {0.80392f,0.40000f,0.00000f}
#define COLORdarkOrange4         {0.54510f,0.27059f,0.00000f}
#define COLORdarkOrchid          {0.60000f,0.19608f,0.80000f}
#define COLORdarkOrchid1         {0.74902f,0.24314f,1.00000f}
#define COLORdarkOrchid2         {0.69804f,0.22745f,0.93333f}
#define COLORdarkOrchid3         {0.60392f,0.19608f,0.80392f}
#define COLORdarkOrchid4         {0.40784f,0.13333f,0.54510f}
#define COLORdarkRed             {0.54510f,0.00000f,0.00000f}
#define COLORdarkSalmon          {0.91373f,0.58824f,0.47843f}
#define COLORdarkSeaGreen        {0.56078f,0.73725f,0.56078f}
#define COLORdarkSeaGreen1       {0.75686f,1.00000f,0.75686f}
#define COLORdarkSeaGreen2       {0.70588f,0.93333f,0.70588f}
#define COLORdarkSeaGreen3       {0.60784f,0.80392f,0.60784f}
#define COLORdarkSeaGreen4       {0.41176f,0.54510f,0.41176f}
#define COLORdarkSlateBlue       {0.28235f,0.23922f,0.54510f}
#define COLORdarkSlateGray       {0.18431f,0.30980f,0.30980f}
#define COLORdarkSlateGray1      {0.59216f,1.00000f,1.00000f}
#define COLORdarkSlateGray2      {0.55294f,0.93333f,0.93333f}
#define COLORdarkSlateGray3      {0.47451f,0.80392f,0.80392f}
#define COLORdarkSlateGray4      {0.32157f,0.54510f,0.54510f}
#define COLORdarkSlateGrey       {0.18431f,0.30980f,0.30980f}
#define COLORdarkTurquoise       {0.00000f,0.80784f,0.81961f}
#define COLORdarkViolet          {0.58039f,0.00000f,0.82745f}
#define COLORdeepPink            {1.00000f,0.07843f,0.57647f}
#define COLORdeepPink1           {1.00000f,0.07843f,0.57647f}
#define COLORdeepPink2           {0.93333f,0.07059f,0.53725f}
#define COLORdeepPink3           {0.80392f,0.06275f,0.46275f}
#define COLORdeepPink4           {0.54510f,0.03922f,0.31373f}
#define COLORdeepSkyBlue         {0.00000f,0.74902f,1.00000f}
#define COLORdeepSkyBlue1        {0.00000f,0.74902f,1.00000f}
#define COLORdeepSkyBlue2        {0.00000f,0.69804f,0.93333f}
#define COLORdeepSkyBlue3        {0.00000f,0.60392f,0.80392f}
#define COLORdeepSkyBlue4        {0.00000f,0.40784f,0.54510f}
#define COLORdimGray             {0.41176f,0.41176f,0.41176f}
#define COLORdimGrey             {0.41176f,0.41176f,0.41176f}
#define COLORdodgerBlue          {0.11765f,0.56471f,1.00000f}
#define COLORdodgerBlue1         {0.11765f,0.56471f,1.00000f}
#define COLORdodgerBlue2         {0.10980f,0.52549f,0.93333f}
#define COLORdodgerBlue3         {0.09412f,0.45490f,0.80392f}
#define COLORdodgerBlue4         {0.06275f,0.30588f,0.54510f}
#define COLORfirebrick           {0.69804f,0.13333f,0.13333f}
#define COLORfirebrick1          {1.00000f,0.18824f,0.18824f}
#define COLORfirebrick2          {0.93333f,0.17255f,0.17255f}
#define COLORfirebrick3          {0.80392f,0.14902f,0.14902f}
#define COLORfirebrick4          {0.54510f,0.10196f,0.10196f}
#define COLORfloralWhite         {1.00000f,0.98039f,0.94118f}
#define COLORforestGreen         {0.13333f,0.54510f,0.13333f}
#define COLORgainsboro           {0.86275f,0.86275f,0.86275f}
#define COLORghostWhite          {0.97255f,0.97255f,1.00000f}
#define COLORgold                {1.00000f,0.84314f,0.00000f}
#define COLORgold1               {1.00000f,0.84314f,0.00000f}
#define COLORgold2               {0.93333f,0.78824f,0.00000f}
#define COLORgold3               {0.80392f,0.67843f,0.00000f}
#define COLORgold4               {0.54510f,0.45882f,0.00000f}
#define COLORgoldenrod           {0.85490f,0.64706f,0.12549f}
#define COLORgoldenrod1          {1.00000f,0.75686f,0.14510f}
#define COLORgoldenrod2          {0.93333f,0.70588f,0.13333f}
#define COLORgoldenrod3          {0.80392f,0.60784f,0.11373f}
#define COLORgoldenrod4          {0.54510f,0.41176f,0.07843f}
#define COLORgray                {0.74510f,0.74510f,0.74510f}
#define COLORgray0               {0.00000f,0.00000f,0.00000f}
#define COLORgray1               {0.01176f,0.01176f,0.01176f}
#define COLORgray10              {0.10196f,0.10196f,0.10196f}
#define COLORgray100             {1.00000f,1.00000f,1.00000f}
#define COLORgray11              {0.10980f,0.10980f,0.10980f}
#define COLORgray12              {0.12157f,0.12157f,0.12157f}
#define COLORgray13              {0.12941f,0.12941f,0.12941f}
#define COLORgray14              {0.14118f,0.14118f,0.14118f}
#define COLORgray15              {0.14902f,0.14902f,0.14902f}
#define COLORgray16              {0.16078f,0.16078f,0.16078f}
#define COLORgray17              {0.16863f,0.16863f,0.16863f}
#define COLORgray18              {0.18039f,0.18039f,0.18039f}
#define COLORgray19              {0.18824f,0.18824f,0.18824f}
#define COLORgray2               {0.01961f,0.01961f,0.01961f}
#define COLORgray20              {0.20000f,0.20000f,0.20000f}
#define COLORgray21              {0.21176f,0.21176f,0.21176f}
#define COLORgray22              {0.21961f,0.21961f,0.21961f}
#define COLORgray23              {0.23137f,0.23137f,0.23137f}
#define COLORgray24              {0.23922f,0.23922f,0.23922f}
#define COLORgray25              {0.25098f,0.25098f,0.25098f}
#define COLORgray26              {0.25882f,0.25882f,0.25882f}
#define COLORgray27              {0.27059f,0.27059f,0.27059f}
#define COLORgray28              {0.27843f,0.27843f,0.27843f}
#define COLORgray29              {0.29020f,0.29020f,0.29020f}
#define COLORgray3               {0.03137f,0.03137f,0.03137f}
#define COLORgray30              {0.30196f,0.30196f,0.30196f}
#define COLORgray31              {0.30980f,0.30980f,0.30980f}
#define COLORgray32              {0.32157f,0.32157f,0.32157f}
#define COLORgray33              {0.32941f,0.32941f,0.32941f}
#define COLORgray34              {0.34118f,0.34118f,0.34118f}
#define COLORgray35              {0.34902f,0.34902f,0.34902f}
#define COLORgray36              {0.36078f,0.36078f,0.36078f}
#define COLORgray37              {0.36863f,0.36863f,0.36863f}
#define COLORgray38              {0.38039f,0.38039f,0.38039f}
#define COLORgray39              {0.38824f,0.38824f,0.38824f}
#define COLORgray4               {0.03922f,0.03922f,0.03922f}
#define COLORgray40              {0.40000f,0.40000f,0.40000f}
#define COLORgray41              {0.41176f,0.41176f,0.41176f}
#define COLORgray42              {0.41961f,0.41961f,0.41961f}
#define COLORgray43              {0.43137f,0.43137f,0.43137f}
#define COLORgray44              {0.43922f,0.43922f,0.43922f}
#define COLORgray45              {0.45098f,0.45098f,0.45098f}
#define COLORgray46              {0.45882f,0.45882f,0.45882f}
#define COLORgray47              {0.47059f,0.47059f,0.47059f}
#define COLORgray48              {0.47843f,0.47843f,0.47843f}
#define COLORgray49              {0.49020f,0.49020f,0.49020f}
#define COLORgray5               {0.05098f,0.05098f,0.05098f}
#define COLORgray50              {0.49804f,0.49804f,0.49804f}
#define COLORgray51              {0.50980f,0.50980f,0.50980f}
#define COLORgray52              {0.52157f,0.52157f,0.52157f}
#define COLORgray53              {0.52941f,0.52941f,0.52941f}
#define COLORgray54              {0.54118f,0.54118f,0.54118f}
#define COLORgray55              {0.54902f,0.54902f,0.54902f}
#define COLORgray56              {0.56078f,0.56078f,0.56078f}
#define COLORgray57              {0.56863f,0.56863f,0.56863f}
#define COLORgray58              {0.58039f,0.58039f,0.58039f}
#define COLORgray59              {0.58824f,0.58824f,0.58824f}
#define COLORgray6               {0.05882f,0.05882f,0.05882f}
#define COLORgray60              {0.60000f,0.60000f,0.60000f}
#define COLORgray61              {0.61176f,0.61176f,0.61176f}
#define COLORgray62              {0.61961f,0.61961f,0.61961f}
#define COLORgray63              {0.63137f,0.63137f,0.63137f}
#define COLORgray64              {0.63922f,0.63922f,0.63922f}
#define COLORgray65              {0.65098f,0.65098f,0.65098f}
#define COLORgray66              {0.65882f,0.65882f,0.65882f}
#define COLORgray67              {0.67059f,0.67059f,0.67059f}
#define COLORgray68              {0.67843f,0.67843f,0.67843f}
#define COLORgray69              {0.69020f,0.69020f,0.69020f}
#define COLORgray7               {0.07059f,0.07059f,0.07059f}
#define COLORgray70              {0.70196f,0.70196f,0.70196f}
#define COLORgray71              {0.70980f,0.70980f,0.70980f}
#define COLORgray72              {0.72157f,0.72157f,0.72157f}
#define COLORgray73              {0.72941f,0.72941f,0.72941f}
#define COLORgray74              {0.74118f,0.74118f,0.74118f}
#define COLORgray75              {0.74902f,0.74902f,0.74902f}
#define COLORgray76              {0.76078f,0.76078f,0.76078f}
#define COLORgray77              {0.76863f,0.76863f,0.76863f}
#define COLORgray78              {0.78039f,0.78039f,0.78039f}
#define COLORgray79              {0.78824f,0.78824f,0.78824f}
#define COLORgray8               {0.07843f,0.07843f,0.07843f}
#define COLORgray80              {0.80000f,0.80000f,0.80000f}
#define COLORgray81              {0.81176f,0.81176f,0.81176f}
#define COLORgray82              {0.81961f,0.81961f,0.81961f}
#define COLORgray83              {0.83137f,0.83137f,0.83137f}
#define COLORgray84              {0.83922f,0.83922f,0.83922f}
#define COLORgray85              {0.85098f,0.85098f,0.85098f}
#define COLORgray86              {0.85882f,0.85882f,0.85882f}
#define COLORgray87              {0.87059f,0.87059f,0.87059f}
#define COLORgray88              {0.87843f,0.87843f,0.87843f}
#define COLORgray89              {0.89020f,0.89020f,0.89020f}
#define COLORgray9               {0.09020f,0.09020f,0.09020f}
#define COLORgray90              {0.89804f,0.89804f,0.89804f}
#define COLORgray91              {0.90980f,0.90980f,0.90980f}
#define COLORgray92              {0.92157f,0.92157f,0.92157f}
#define COLORgray93              {0.92941f,0.92941f,0.92941f}
#define COLORgray94              {0.94118f,0.94118f,0.94118f}
#define COLORgray95              {0.94902f,0.94902f,0.94902f}
#define COLORgray96              {0.96078f,0.96078f,0.96078f}
#define COLORgray97              {0.96863f,0.96863f,0.96863f}
#define COLORgray98              {0.98039f,0.98039f,0.98039f}
#define COLORgray99              {0.98824f,0.98824f,0.98824f}
#define COLORgreen               {0.00000f,1.00000f,0.00000f}
#define COLORgreen1              {0.00000f,1.00000f,0.00000f}
#define COLORgreen2              {0.00000f,0.93333f,0.00000f}
#define COLORgreen3              {0.00000f,0.80392f,0.00000f}
#define COLORgreen4              {0.00000f,0.54510f,0.00000f}
#define COLORgreenYellow         {0.67843f,1.00000f,0.18431f}
#define COLORgrey                {0.74510f,0.74510f,0.74510f}
#define COLORgrey0               {0.00000f,0.00000f,0.00000f}
#define COLORgrey1               {0.01176f,0.01176f,0.01176f}
#define COLORgrey10              {0.10196f,0.10196f,0.10196f}
#define COLORgrey100             {1.00000f,1.00000f,1.00000f}
#define COLORgrey11              {0.10980f,0.10980f,0.10980f}
#define COLORgrey12              {0.12157f,0.12157f,0.12157f}
#define COLORgrey13              {0.12941f,0.12941f,0.12941f}
#define COLORgrey14              {0.14118f,0.14118f,0.14118f}
#define COLORgrey15              {0.14902f,0.14902f,0.14902f}
#define COLORgrey16              {0.16078f,0.16078f,0.16078f}
#define COLORgrey17              {0.16863f,0.16863f,0.16863f}
#define COLORgrey18              {0.18039f,0.18039f,0.18039f}
#define COLORgrey19              {0.18824f,0.18824f,0.18824f}
#define COLORgrey2               {0.01961f,0.01961f,0.01961f}
#define COLORgrey20              {0.20000f,0.20000f,0.20000f}
#define COLORgrey21              {0.21176f,0.21176f,0.21176f}
#define COLORgrey22              {0.21961f,0.21961f,0.21961f}
#define COLORgrey23              {0.23137f,0.23137f,0.23137f}
#define COLORgrey24              {0.23922f,0.23922f,0.23922f}
#define COLORgrey25              {0.25098f,0.25098f,0.25098f}
#define COLORgrey26              {0.25882f,0.25882f,0.25882f}
#define COLORgrey27              {0.27059f,0.27059f,0.27059f}
#define COLORgrey28              {0.27843f,0.27843f,0.27843f}
#define COLORgrey29              {0.29020f,0.29020f,0.29020f}
#define COLORgrey3               {0.03137f,0.03137f,0.03137f}
#define COLORgrey30              {0.30196f,0.30196f,0.30196f}
#define COLORgrey31              {0.30980f,0.30980f,0.30980f}
#define COLORgrey32              {0.32157f,0.32157f,0.32157f}
#define COLORgrey33              {0.32941f,0.32941f,0.32941f}
#define COLORgrey34              {0.34118f,0.34118f,0.34118f}
#define COLORgrey35              {0.34902f,0.34902f,0.34902f}
#define COLORgrey36              {0.36078f,0.36078f,0.36078f}
#define COLORgrey37              {0.36863f,0.36863f,0.36863f}
#define COLORgrey38              {0.38039f,0.38039f,0.38039f}
#define COLORgrey39              {0.38824f,0.38824f,0.38824f}
#define COLORgrey4               {0.03922f,0.03922f,0.03922f}
#define COLORgrey40              {0.40000f,0.40000f,0.40000f}
#define COLORgrey41              {0.41176f,0.41176f,0.41176f}
#define COLORgrey42              {0.41961f,0.41961f,0.41961f}
#define COLORgrey43              {0.43137f,0.43137f,0.43137f}
#define COLORgrey44              {0.43922f,0.43922f,0.43922f}
#define COLORgrey45              {0.45098f,0.45098f,0.45098f}
#define COLORgrey46              {0.45882f,0.45882f,0.45882f}
#define COLORgrey47              {0.47059f,0.47059f,0.47059f}
#define COLORgrey48              {0.47843f,0.47843f,0.47843f}
#define COLORgrey49              {0.49020f,0.49020f,0.49020f}
#define COLORgrey5               {0.05098f,0.05098f,0.05098f}
#define COLORgrey50              {0.49804f,0.49804f,0.49804f}
#define COLORgrey51              {0.50980f,0.50980f,0.50980f}
#define COLORgrey52              {0.52157f,0.52157f,0.52157f}
#define COLORgrey53              {0.52941f,0.52941f,0.52941f}
#define COLORgrey54              {0.54118f,0.54118f,0.54118f}
#define COLORgrey55              {0.54902f,0.54902f,0.54902f}
#define COLORgrey56              {0.56078f,0.56078f,0.56078f}
#define COLORgrey57              {0.56863f,0.56863f,0.56863f}
#define COLORgrey58              {0.58039f,0.58039f,0.58039f}
#define COLORgrey59              {0.58824f,0.58824f,0.58824f}
#define COLORgrey6               {0.05882f,0.05882f,0.05882f}
#define COLORgrey60              {0.60000f,0.60000f,0.60000f}
#define COLORgrey61              {0.61176f,0.61176f,0.61176f}
#define COLORgrey62              {0.61961f,0.61961f,0.61961f}
#define COLORgrey63              {0.63137f,0.63137f,0.63137f}
#define COLORgrey64              {0.63922f,0.63922f,0.63922f}
#define COLORgrey65              {0.65098f,0.65098f,0.65098f}
#define COLORgrey66              {0.65882f,0.65882f,0.65882f}
#define COLORgrey67              {0.67059f,0.67059f,0.67059f}
#define COLORgrey68              {0.67843f,0.67843f,0.67843f}
#define COLORgrey69              {0.69020f,0.69020f,0.69020f}
#define COLORgrey7               {0.07059f,0.07059f,0.07059f}
#define COLORgrey70              {0.70196f,0.70196f,0.70196f}
#define COLORgrey71              {0.70980f,0.70980f,0.70980f}
#define COLORgrey72              {0.72157f,0.72157f,0.72157f}
#define COLORgrey73              {0.72941f,0.72941f,0.72941f}
#define COLORgrey74              {0.74118f,0.74118f,0.74118f}
#define COLORgrey75              {0.74902f,0.74902f,0.74902f}
#define COLORgrey76              {0.76078f,0.76078f,0.76078f}
#define COLORgrey77              {0.76863f,0.76863f,0.76863f}
#define COLORgrey78              {0.78039f,0.78039f,0.78039f}
#define COLORgrey79              {0.78824f,0.78824f,0.78824f}
#define COLORgrey8               {0.07843f,0.07843f,0.07843f}
#define COLORgrey80              {0.80000f,0.80000f,0.80000f}
#define COLORgrey81              {0.81176f,0.81176f,0.81176f}
#define COLORgrey82              {0.81961f,0.81961f,0.81961f}
#define COLORgrey83              {0.83137f,0.83137f,0.83137f}
#define COLORgrey84              {0.83922f,0.83922f,0.83922f}
#define COLORgrey85              {0.85098f,0.85098f,0.85098f}
#define COLORgrey86              {0.85882f,0.85882f,0.85882f}
#define COLORgrey87              {0.87059f,0.87059f,0.87059f}
#define COLORgrey88              {0.87843f,0.87843f,0.87843f}
#define COLORgrey89              {0.89020f,0.89020f,0.89020f}
#define COLORgrey9               {0.09020f,0.09020f,0.09020f}
#define COLORgrey90              {0.89804f,0.89804f,0.89804f}
#define COLORgrey91              {0.90980f,0.90980f,0.90980f}
#define COLORgrey92              {0.92157f,0.92157f,0.92157f}
#define COLORgrey93              {0.92941f,0.92941f,0.92941f}
#define COLORgrey94              {0.94118f,0.94118f,0.94118f}
#define COLORgrey95              {0.94902f,0.94902f,0.94902f}
#define COLORgrey96              {0.96078f,0.96078f,0.96078f}
#define COLORgrey97              {0.96863f,0.96863f,0.96863f}
#define COLORgrey98              {0.98039f,0.98039f,0.98039f}
#define COLORgrey99              {0.98824f,0.98824f,0.98824f}
#define COLORhoneydew            {0.94118f,1.00000f,0.94118f}
#define COLORhoneydew1           {0.94118f,1.00000f,0.94118f}
#define COLORhoneydew2           {0.87843f,0.93333f,0.87843f}
#define COLORhoneydew3           {0.75686f,0.80392f,0.75686f}
#define COLORhoneydew4           {0.51373f,0.54510f,0.51373f}
#define COLORhotPink             {1.00000f,0.41176f,0.70588f}
#define COLORhotPink1            {1.00000f,0.43137f,0.70588f}
#define COLORhotPink2            {0.93333f,0.41569f,0.65490f}
#define COLORhotPink3            {0.80392f,0.37647f,0.56471f}
#define COLORhotPink4            {0.54510f,0.22745f,0.38431f}
#define COLORindianRed           {0.80392f,0.36078f,0.36078f}
#define COLORindianRed1          {1.00000f,0.41569f,0.41569f}
#define COLORindianRed2          {0.93333f,0.38824f,0.38824f}
#define COLORindianRed3          {0.80392f,0.33333f,0.33333f}
#define COLORindianRed4          {0.54510f,0.22745f,0.22745f}
#define COLORivory               {1.00000f,1.00000f,0.94118f}
#define COLORivory1              {1.00000f,1.00000f,0.94118f}
#define COLORivory2              {0.93333f,0.93333f,0.87843f}
#define COLORivory3              {0.80392f,0.80392f,0.75686f}
#define COLORivory4              {0.54510f,0.54510f,0.51373f}
#define COLORkhaki               {0.94118f,0.90196f,0.54902f}
#define COLORkhaki1              {1.00000f,0.96471f,0.56078f}
#define COLORkhaki2              {0.93333f,0.90196f,0.52157f}
#define COLORkhaki3              {0.80392f,0.77647f,0.45098f}
#define COLORkhaki4              {0.54510f,0.52549f,0.30588f}
#define COLORlavender            {0.90196f,0.90196f,0.98039f}
#define COLORlavenderBlush       {1.00000f,0.94118f,0.96078f}
#define COLORlavenderBlush1      {1.00000f,0.94118f,0.96078f}
#define COLORlavenderBlush2      {0.93333f,0.87843f,0.89804f}
#define COLORlavenderBlush3      {0.80392f,0.75686f,0.77255f}
#define COLORlavenderBlush4      {0.54510f,0.51373f,0.52549f}
#define COLORlawnGreen           {0.48627f,0.98824f,0.00000f}
#define COLORlemonChiffon        {1.00000f,0.98039f,0.80392f}
#define COLORlemonChiffon1       {1.00000f,0.98039f,0.80392f}
#define COLORlemonChiffon2       {0.93333f,0.91373f,0.74902f}
#define COLORlemonChiffon3       {0.80392f,0.78824f,0.64706f}
#define COLORlemonChiffon4       {0.54510f,0.53725f,0.43922f}
#define COLORlightBlue           {0.67843f,0.84706f,0.90196f}
#define COLORlightBlue1          {0.74902f,0.93725f,1.00000f}
#define COLORlightBlue2          {0.69804f,0.87451f,0.93333f}
#define COLORlightBlue3          {0.60392f,0.75294f,0.80392f}
#define COLORlightBlue4          {0.40784f,0.51373f,0.54510f}
#define COLORlightCoral          {0.94118f,0.50196f,0.50196f}
#define COLORlightCyan           {0.87843f,1.00000f,1.00000f}
#define COLORlightCyan1          {0.87843f,1.00000f,1.00000f}
#define COLORlightCyan2          {0.81961f,0.93333f,0.93333f}
#define COLORlightCyan3          {0.70588f,0.80392f,0.80392f}
#define COLORlightCyan4          {0.47843f,0.54510f,0.54510f}
#define COLORlightGoldenrod      {0.93333f,0.86667f,0.50980f}
#define COLORlightGoldenrod1     {1.00000f,0.92549f,0.54510f}
#define COLORlightGoldenrod2     {0.93333f,0.86275f,0.50980f}
#define COLORlightGoldenrod3     {0.80392f,0.74510f,0.43922f}
#define COLORlightGoldenrod4     {0.54510f,0.50588f,0.29804f}
#define COLORlightGoldenrodYellow {0.98039f,0.98039f,0.82353f}
#define COLORlightGray           {0.82745f,0.82745f,0.82745f}
#define COLORlightGreen          {0.56471f,0.93333f,0.56471f}
#define COLORlightGrey           {0.82745f,0.82745f,0.82745f}
#define COLORlightPink           {1.00000f,0.71373f,0.75686f}
#define COLORlightPink1          {1.00000f,0.68235f,0.72549f}
#define COLORlightPink2          {0.93333f,0.63529f,0.67843f}
#define COLORlightPink3          {0.80392f,0.54902f,0.58431f}
#define COLORlightPink4          {0.54510f,0.37255f,0.39608f}
#define COLORlightSalmon         {1.00000f,0.62745f,0.47843f}
#define COLORlightSalmon1        {1.00000f,0.62745f,0.47843f}
#define COLORlightSalmon2        {0.93333f,0.58431f,0.44706f}
#define COLORlightSalmon3        {0.80392f,0.50588f,0.38431f}
#define COLORlightSalmon4        {0.54510f,0.34118f,0.25882f}
#define COLORlightSeaGreen       {0.12549f,0.69804f,0.66667f}
#define COLORlightSkyBlue        {0.52941f,0.80784f,0.98039f}
#define COLORlightSkyBlue1       {0.69020f,0.88627f,1.00000f}
#define COLORlightSkyBlue2       {0.64314f,0.82745f,0.93333f}
#define COLORlightSkyBlue3       {0.55294f,0.71373f,0.80392f}
#define COLORlightSkyBlue4       {0.37647f,0.48235f,0.54510f}
#define COLORlightSlateBlue      {0.51765f,0.43922f,1.00000f}
#define COLORlightSlateGray      {0.46667f,0.53333f,0.60000f}
#define COLORlightSlateGrey      {0.46667f,0.53333f,0.60000f}
#define COLORlightSteelBlue      {0.69020f,0.76863f,0.87059f}
#define COLORlightSteelBlue1     {0.79216f,0.88235f,1.00000f}
#define COLORlightSteelBlue2     {0.73725f,0.82353f,0.93333f}
#define COLORlightSteelBlue3     {0.63529f,0.70980f,0.80392f}
#define COLORlightSteelBlue4     {0.43137f,0.48235f,0.54510f}
#define COLORlightYellow         {1.00000f,1.00000f,0.87843f}
#define COLORlightYellow1        {1.00000f,1.00000f,0.87843f}
#define COLORlightYellow2        {0.93333f,0.93333f,0.81961f}
#define COLORlightYellow3        {0.80392f,0.80392f,0.70588f}
#define COLORlightYellow4        {0.54510f,0.54510f,0.47843f}
#define COLORlimeGreen           {0.19608f,0.80392f,0.19608f}
#define COLORlinen               {0.98039f,0.94118f,0.90196f}
#define COLORmagenta             {1.00000f,0.00000f,1.00000f}
#define COLORmagenta1            {1.00000f,0.00000f,1.00000f}
#define COLORmagenta2            {0.93333f,0.00000f,0.93333f}
#define COLORmagenta3            {0.80392f,0.00000f,0.80392f}
#define COLORmagenta4            {0.54510f,0.00000f,0.54510f}
#define COLORmaroon              {0.69020f,0.18824f,0.37647f}
#define COLORmaroon1             {1.00000f,0.20392f,0.70196f}
#define COLORmaroon2             {0.93333f,0.18824f,0.65490f}
#define COLORmaroon3             {0.80392f,0.16078f,0.56471f}
#define COLORmaroon4             {0.54510f,0.10980f,0.38431f}
#define COLORmediumAquamarine    {0.40000f,0.80392f,0.66667f}
#define COLORmediumBlue          {0.00000f,0.00000f,0.80392f}
#define COLORmediumOrchid        {0.72941f,0.33333f,0.82745f}
#define COLORmediumOrchid1       {0.87843f,0.40000f,1.00000f}
#define COLORmediumOrchid2       {0.81961f,0.37255f,0.93333f}
#define COLORmediumOrchid3       {0.70588f,0.32157f,0.80392f}
#define COLORmediumOrchid4       {0.47843f,0.21569f,0.54510f}
#define COLORmediumPurple        {0.57647f,0.43922f,0.85882f}
#define COLORmediumPurple1       {0.67059f,0.50980f,1.00000f}
#define COLORmediumPurple2       {0.62353f,0.47451f,0.93333f}
#define COLORmediumPurple3       {0.53725f,0.40784f,0.80392f}
#define COLORmediumPurple4       {0.36471f,0.27843f,0.54510f}
#define COLORmediumSeaGreen      {0.23529f,0.70196f,0.44314f}
#define COLORmediumSlateBlue     {0.48235f,0.40784f,0.93333f}
#define COLORmediumSpringGreen   {0.00000f,0.98039f,0.60392f}
#define COLORmediumTurquoise     {0.28235f,0.81961f,0.80000f}
#define COLORmediumVioletRed     {0.78039f,0.08235f,0.52157f}
#define COLORmidnightBlue        {0.09804f,0.09804f,0.43922f}
#define COLORmintCream           {0.96078f,1.00000f,0.98039f}
#define COLORmistyRose           {1.00000f,0.89412f,0.88235f}
#define COLORmistyRose1          {1.00000f,0.89412f,0.88235f}
#define COLORmistyRose2          {0.93333f,0.83529f,0.82353f}
#define COLORmistyRose3          {0.80392f,0.71765f,0.70980f}
#define COLORmistyRose4          {0.54510f,0.49020f,0.48235f}
#define COLORmoccasin            {1.00000f,0.89412f,0.70980f}
#define COLORnavajoWhite         {1.00000f,0.87059f,0.67843f}
#define COLORnavajoWhite1        {1.00000f,0.87059f,0.67843f}
#define COLORnavajoWhite2        {0.93333f,0.81176f,0.63137f}
#define COLORnavajoWhite3        {0.80392f,0.70196f,0.54510f}
#define COLORnavajoWhite4        {0.54510f,0.47451f,0.36863f}
#define COLORnavy                {0.00000f,0.00000f,0.50196f}
#define COLORnavyBlue            {0.00000f,0.00000f,0.50196f}
#define COLORoldLace             {0.99216f,0.96078f,0.90196f}
#define COLORoliveDrab           {0.41961f,0.55686f,0.13725f}
#define COLORoliveDrab1          {0.75294f,1.00000f,0.24314f}
#define COLORoliveDrab2          {0.70196f,0.93333f,0.22745f}
#define COLORoliveDrab3          {0.60392f,0.80392f,0.19608f}
#define COLORoliveDrab4          {0.41176f,0.54510f,0.13333f}
#define COLORorange              {1.00000f,0.64706f,0.00000f}
#define COLORorange1             {1.00000f,0.64706f,0.00000f}
#define COLORorange2             {0.93333f,0.60392f,0.00000f}
#define COLORorange3             {0.80392f,0.52157f,0.00000f}
#define COLORorange4             {0.54510f,0.35294f,0.00000f}
#define COLORorangeRed           {1.00000f,0.27059f,0.00000f}
#define COLORorangeRed1          {1.00000f,0.27059f,0.00000f}
#define COLORorangeRed2          {0.93333f,0.25098f,0.00000f}
#define COLORorangeRed3          {0.80392f,0.21569f,0.00000f}
#define COLORorangeRed4          {0.54510f,0.14510f,0.00000f}
#define COLORorchid              {0.85490f,0.43922f,0.83922f}
#define COLORorchid1             {1.00000f,0.51373f,0.98039f}
#define COLORorchid2             {0.93333f,0.47843f,0.91373f}
#define COLORorchid3             {0.80392f,0.41176f,0.78824f}
#define COLORorchid4             {0.54510f,0.27843f,0.53725f}
#define COLORpaleGoldenrod       {0.93333f,0.90980f,0.66667f}
#define COLORpaleGreen           {0.59608f,0.98431f,0.59608f}
#define COLORpaleGreen1          {0.60392f,1.00000f,0.60392f}
#define COLORpaleGreen2          {0.56471f,0.93333f,0.56471f}
#define COLORpaleGreen3          {0.48627f,0.80392f,0.48627f}
#define COLORpaleGreen4          {0.32941f,0.54510f,0.32941f}
#define COLORpaleTurquoise       {0.68627f,0.93333f,0.93333f}
#define COLORpaleTurquoise1      {0.73333f,1.00000f,1.00000f}
#define COLORpaleTurquoise2      {0.68235f,0.93333f,0.93333f}
#define COLORpaleTurquoise3      {0.58824f,0.80392f,0.80392f}
#define COLORpaleTurquoise4      {0.40000f,0.54510f,0.54510f}
#define COLORpaleVioletRed       {0.85882f,0.43922f,0.57647f}
#define COLORpaleVioletRed1      {1.00000f,0.50980f,0.67059f}
#define COLORpaleVioletRed2      {0.93333f,0.47451f,0.62353f}
#define COLORpaleVioletRed3      {0.80392f,0.40784f,0.53725f}
#define COLORpaleVioletRed4      {0.54510f,0.27843f,0.36471f}
#define COLORpapayaWhip          {1.00000f,0.93725f,0.83529f}
#define COLORpeachPuff           {1.00000f,0.85490f,0.72549f}
#define COLORpeachPuff1          {1.00000f,0.85490f,0.72549f}
#define COLORpeachPuff2          {0.93333f,0.79608f,0.67843f}
#define COLORpeachPuff3          {0.80392f,0.68627f,0.58431f}
#define COLORpeachPuff4          {0.54510f,0.46667f,0.39608f}
#define COLORperu                {0.80392f,0.52157f,0.24706f}
#define COLORpink                {1.00000f,0.75294f,0.79608f}
#define COLORpink1               {1.00000f,0.70980f,0.77255f}
#define COLORpink2               {0.93333f,0.66275f,0.72157f}
#define COLORpink3               {0.80392f,0.56863f,0.61961f}
#define COLORpink4               {0.54510f,0.38824f,0.42353f}
#define COLORplum                {0.86667f,0.62745f,0.86667f}
#define COLORplum1               {1.00000f,0.73333f,1.00000f}
#define COLORplum2               {0.93333f,0.68235f,0.93333f}
#define COLORplum3               {0.80392f,0.58824f,0.80392f}
#define COLORplum4               {0.54510f,0.40000f,0.54510f}
#define COLORpowderBlue          {0.69020f,0.87843f,0.90196f}
#define COLORpurple              {0.62745f,0.12549f,0.94118f}
#define COLORpurple1             {0.60784f,0.18824f,1.00000f}
#define COLORpurple2             {0.56863f,0.17255f,0.93333f}
#define COLORpurple3             {0.49020f,0.14902f,0.80392f}
#define COLORpurple4             {0.33333f,0.10196f,0.54510f}
#define COLORred                 {1.00000f,0.00000f,0.00000f}
#define COLORred1                {1.00000f,0.00000f,0.00000f}
#define COLORred2                {0.93333f,0.00000f,0.00000f}
#define COLORred3                {0.80392f,0.00000f,0.00000f}
#define COLORred4                {0.54510f,0.00000f,0.00000f}
#define COLORrosyBrown           {0.73725f,0.56078f,0.56078f}
#define COLORrosyBrown1          {1.00000f,0.75686f,0.75686f}
#define COLORrosyBrown2          {0.93333f,0.70588f,0.70588f}
#define COLORrosyBrown3          {0.80392f,0.60784f,0.60784f}
#define COLORrosyBrown4          {0.54510f,0.41176f,0.41176f}
#define COLORroyalBlue           {0.25490f,0.41176f,0.88235f}
#define COLORroyalBlue1          {0.28235f,0.46275f,1.00000f}
#define COLORroyalBlue2          {0.26275f,0.43137f,0.93333f}
#define COLORroyalBlue3          {0.22745f,0.37255f,0.80392f}
#define COLORroyalBlue4          {0.15294f,0.25098f,0.54510f}
#define COLORsaddleBrown         {0.54510f,0.27059f,0.07451f}
#define COLORsalmon              {0.98039f,0.50196f,0.44706f}
#define COLORsalmon1             {1.00000f,0.54902f,0.41176f}
#define COLORsalmon2             {0.93333f,0.50980f,0.38431f}
#define COLORsalmon3             {0.80392f,0.43922f,0.32941f}
#define COLORsalmon4             {0.54510f,0.29804f,0.22353f}
#define COLORsandyBrown          {0.95686f,0.64314f,0.37647f}
#define COLORseaGreen            {0.18039f,0.54510f,0.34118f}
#define COLORseaGreen1           {0.32941f,1.00000f,0.62353f}
#define COLORseaGreen2           {0.30588f,0.93333f,0.58039f}
#define COLORseaGreen3           {0.26275f,0.80392f,0.50196f}
#define COLORseaGreen4           {0.18039f,0.54510f,0.34118f}
#define COLORseashell            {1.00000f,0.96078f,0.93333f}
#define COLORseashell1           {1.00000f,0.96078f,0.93333f}
#define COLORseashell2           {0.93333f,0.89804f,0.87059f}
#define COLORseashell3           {0.80392f,0.77255f,0.74902f}
#define COLORseashell4           {0.54510f,0.52549f,0.50980f}
#define COLORsienna              {0.62745f,0.32157f,0.17647f}
#define COLORsienna1             {1.00000f,0.50980f,0.27843f}
#define COLORsienna2             {0.93333f,0.47451f,0.25882f}
#define COLORsienna3             {0.80392f,0.40784f,0.22353f}
#define COLORsienna4             {0.54510f,0.27843f,0.14902f}
#define COLORskyBlue             {0.52941f,0.80784f,0.92157f}
#define COLORskyBlue1            {0.52941f,0.80784f,1.00000f}
#define COLORskyBlue2            {0.49412f,0.75294f,0.93333f}
#define COLORskyBlue3            {0.42353f,0.65098f,0.80392f}
#define COLORskyBlue4            {0.29020f,0.43922f,0.54510f}
#define COLORslateBlue           {0.41569f,0.35294f,0.80392f}
#define COLORslateBlue1          {0.51373f,0.43529f,1.00000f}
#define COLORslateBlue2          {0.47843f,0.40392f,0.93333f}
#define COLORslateBlue3          {0.41176f,0.34902f,0.80392f}
#define COLORslateBlue4          {0.27843f,0.23529f,0.54510f}
#define COLORslateGray           {0.43922f,0.50196f,0.56471f}
#define COLORslateGray1          {0.77647f,0.88627f,1.00000f}
#define COLORslateGray2          {0.72549f,0.82745f,0.93333f}
#define COLORslateGray3          {0.62353f,0.71373f,0.80392f}
#define COLORslateGray4          {0.42353f,0.48235f,0.54510f}
#define COLORslateGrey           {0.43922f,0.50196f,0.56471f}
#define COLORsnow                {1.00000f,0.98039f,0.98039f}
#define COLORsnow1               {1.00000f,0.98039f,0.98039f}
#define COLORsnow2               {0.93333f,0.91373f,0.91373f}
#define COLORsnow3               {0.80392f,0.78824f,0.78824f}
#define COLORsnow4               {0.54510f,0.53725f,0.53725f}
#define COLORspringGreen         {0.00000f,1.00000f,0.49804f}
#define COLORspringGreen1        {0.00000f,1.00000f,0.49804f}
#define COLORspringGreen2        {0.00000f,0.93333f,0.46275f}
#define COLORspringGreen3        {0.00000f,0.80392f,0.40000f}
#define COLORspringGreen4        {0.00000f,0.54510f,0.27059f}
#define COLORsteelBlue           {0.27451f,0.50980f,0.70588f}
#define COLORsteelBlue1          {0.38824f,0.72157f,1.00000f}
#define COLORsteelBlue2          {0.36078f,0.67451f,0.93333f}
#define COLORsteelBlue3          {0.30980f,0.58039f,0.80392f}
#define COLORsteelBlue4          {0.21176f,0.39216f,0.54510f}
#define COLORtan1                {1.00000f,0.64706f,0.30980f}
#define COLORtan2                {0.93333f,0.60392f,0.28627f}
#define COLORtan3                {0.80392f,0.52157f,0.24706f}
#define COLORtan4                {0.54510f,0.35294f,0.16863f}
#define COLORthistle             {0.84706f,0.74902f,0.84706f}
#define COLORthistle1            {1.00000f,0.88235f,1.00000f}
#define COLORthistle2            {0.93333f,0.82353f,0.93333f}
#define COLORthistle3            {0.80392f,0.70980f,0.80392f}
#define COLORthistle4            {0.54510f,0.48235f,0.54510f}
#define COLORtomato              {1.00000f,0.38824f,0.27843f}
#define COLORtomato1             {1.00000f,0.38824f,0.27843f}
#define COLORtomato2             {0.93333f,0.36078f,0.25882f}
#define COLORtomato3             {0.80392f,0.30980f,0.22353f}
#define COLORtomato4             {0.54510f,0.21176f,0.14902f}
#define COLORturquoise           {0.25098f,0.87843f,0.81569f}
#define COLORturquoise1          {0.00000f,0.96078f,1.00000f}
#define COLORturquoise2          {0.00000f,0.89804f,0.93333f}
#define COLORturquoise3          {0.00000f,0.77255f,0.80392f}
#define COLORturquoise4          {0.00000f,0.52549f,0.54510f}
#define COLORviolet              {0.93333f,0.50980f,0.93333f}
#define COLORvioletRed           {0.81569f,0.12549f,0.56471f}
#define COLORvioletRed1          {1.00000f,0.24314f,0.58824f}
#define COLORvioletRed2          {0.93333f,0.22745f,0.54902f}
#define COLORvioletRed3          {0.80392f,0.19608f,0.47059f}
#define COLORvioletRed4          {0.54510f,0.13333f,0.32157f}
#define COLORwheat               {0.96078f,0.87059f,0.70196f}
#define COLORwheat1              {1.00000f,0.90588f,0.72941f}
#define COLORwheat2              {0.93333f,0.84706f,0.68235f}
#define COLORwheat3              {0.80392f,0.72941f,0.58824f}
#define COLORwheat4              {0.54510f,0.49412f,0.40000f}
#define COLORwhite               {1.00000f,1.00000f,1.00000f}
#define COLORwhiteSmoke          {0.96078f,0.96078f,0.96078f}
#define COLORyellow              {1.00000f,1.00000f,0.00000f}
#define COLORyellow1             {1.00000f,1.00000f,0.00000f}
#define COLORyellow2             {0.93333f,0.93333f,0.00000f}
#define COLORyellow3             {0.80392f,0.80392f,0.00000f}
#define COLORyellow4             {0.54510f,0.54510f,0.00000f}
#define COLORyellowGreen         {0.60392f,0.80392f,0.19608f}

char *rgbName[656] =
{
	"aliceBlue",
	"antiqueWhite",
	"antiqueWhite1",
	"antiqueWhite2",
	"antiqueWhite3",
	"antiqueWhite4",
	"aquamarine",
	"aquamarine1",
	"aquamarine2",
	"aquamarine3",
	"aquamarine4",
	"azure",
	"azure1",
	"azure2",
	"azure3",
	"azure4",
	"beige",
	"bisque",
	"bisque1",
	"bisque2",
	"bisque3",
	"bisque4",
	"black",
	"blanchedAlmond",
	"blue",
	"blue1",
	"blue2",
	"blue3",
	"blue4",
	"blueViolet",
	"brown",
	"brown1",
	"brown2",
	"brown3",
	"brown4",
	"burlywood",
	"burlywood1",
	"burlywood2",
	"burlywood3",
	"burlywood4",
	"cadetBlue",
	"cadetBlue1",
	"cadetBlue2",
	"cadetBlue3",
	"cadetBlue4",
	"chartreuse",
	"chartreuse1",
	"chartreuse2",
	"chartreuse3",
	"chartreuse4",
	"chocolate",
	"chocolate1",
	"chocolate2",
	"chocolate3",
	"chocolate4",
	"coral",
	"coral1",
	"coral2",
	"coral3",
	"coral4",
	"cornflowerBlue",
	"cornsilk",
	"cornsilk1",
	"cornsilk2",
	"cornsilk3",
	"cornsilk4",
	"cyan",
	"cyan1",
	"cyan2",
	"cyan3",
	"cyan4",
	"darkBlue",
	"darkCyan",
	"darkGoldenrod",
	"darkGoldenrod1",
	"darkGoldenrod2",
	"darkGoldenrod3",
	"darkGoldenrod4",
	"darkGray",
	"darkGreen",
	"darkGrey",
	"darkKhaki",
	"darkMagenta",
	"darkOliveGreen",
	"darkOliveGreen1",
	"darkOliveGreen2",
	"darkOliveGreen3",
	"darkOliveGreen4",
	"darkOrange",
	"darkOrange1",
	"darkOrange2",
	"darkOrange3",
	"darkOrange4",
	"darkOrchid",
	"darkOrchid1",
	"darkOrchid2",
	"darkOrchid3",
	"darkOrchid4",
	"darkRed",
	"darkSalmon",
	"darkSeaGreen",
	"darkSeaGreen1",
	"darkSeaGreen2",
	"darkSeaGreen3",
	"darkSeaGreen4",
	"darkSlateBlue",
	"darkSlateGray",
	"darkSlateGray1",
	"darkSlateGray2",
	"darkSlateGray3",
	"darkSlateGray4",
	"darkSlateGrey",
	"darkTurquoise",
	"darkViolet",
	"deepPink",
	"deepPink1",
	"deepPink2",
	"deepPink3",
	"deepPink4",
	"deepSkyBlue",
	"deepSkyBlue1",
	"deepSkyBlue2",
	"deepSkyBlue3",
	"deepSkyBlue4",
	"dimGray",
	"dimGrey",
	"dodgerBlue",
	"dodgerBlue1",
	"dodgerBlue2",
	"dodgerBlue3",
	"dodgerBlue4",
	"firebrick",
	"firebrick1",
	"firebrick2",
	"firebrick3",
	"firebrick4",
	"floralWhite",
	"forestGreen",
	"gainsboro",
	"ghostWhite",
	"gold",
	"gold1",
	"gold2",
	"gold3",
	"gold4",
	"goldenrod",
	"goldenrod1",
	"goldenrod2",
	"goldenrod3",
	"goldenrod4",
	"gray",
	"gray0",
	"gray1",
	"gray10",
	"gray100",
	"gray11",
	"gray12",
	"gray13",
	"gray14",
	"gray15",
	"gray16",
	"gray17",
	"gray18",
	"gray19",
	"gray2",
	"gray20",
	"gray21",
	"gray22",
	"gray23",
	"gray24",
	"gray25",
	"gray26",
	"gray27",
	"gray28",
	"gray29",
	"gray3",
	"gray30",
	"gray31",
	"gray32",
	"gray33",
	"gray34",
	"gray35",
	"gray36",
	"gray37",
	"gray38",
	"gray39",
	"gray4",
	"gray40",
	"gray41",
	"gray42",
	"gray43",
	"gray44",
	"gray45",
	"gray46",
	"gray47",
	"gray48",
	"gray49",
	"gray5",
	"gray50",
	"gray51",
	"gray52",
	"gray53",
	"gray54",
	"gray55",
	"gray56",
	"gray57",
	"gray58",
	"gray59",
	"gray6",
	"gray60",
	"gray61",
	"gray62",
	"gray63",
	"gray64",
	"gray65",
	"gray66",
	"gray67",
	"gray68",
	"gray69",
	"gray7",
	"gray70",
	"gray71",
	"gray72",
	"gray73",
	"gray74",
	"gray75",
	"gray76",
	"gray77",
	"gray78",
	"gray79",
	"gray8",
	"gray80",
	"gray81",
	"gray82",
	"gray83",
	"gray84",
	"gray85",
	"gray86",
	"gray87",
	"gray88",
	"gray89",
	"gray9",
	"gray90",
	"gray91",
	"gray92",
	"gray93",
	"gray94",
	"gray95",
	"gray96",
	"gray97",
	"gray98",
	"gray99",
	"green",
	"green1",
	"green2",
	"green3",
	"green4",
	"greenYellow",
	"grey",
	"grey0",
	"grey1",
	"grey10",
	"grey100",
	"grey11",
	"grey12",
	"grey13",
	"grey14",
	"grey15",
	"grey16",
	"grey17",
	"grey18",
	"grey19",
	"grey2",
	"grey20",
	"grey21",
	"grey22",
	"grey23",
	"grey24",
	"grey25",
	"grey26",
	"grey27",
	"grey28",
	"grey29",
	"grey3",
	"grey30",
	"grey31",
	"grey32",
	"grey33",
	"grey34",
	"grey35",
	"grey36",
	"grey37",
	"grey38",
	"grey39",
	"grey4",
	"grey40",
	"grey41",
	"grey42",
	"grey43",
	"grey44",
	"grey45",
	"grey46",
	"grey47",
	"grey48",
	"grey49",
	"grey5",
	"grey50",
	"grey51",
	"grey52",
	"grey53",
	"grey54",
	"grey55",
	"grey56",
	"grey57",
	"grey58",
	"grey59",
	"grey6",
	"grey60",
	"grey61",
	"grey62",
	"grey63",
	"grey64",
	"grey65",
	"grey66",
	"grey67",
	"grey68",
	"grey69",
	"grey7",
	"grey70",
	"grey71",
	"grey72",
	"grey73",
	"grey74",
	"grey75",
	"grey76",
	"grey77",
	"grey78",
	"grey79",
	"grey8",
	"grey80",
	"grey81",
	"grey82",
	"grey83",
	"grey84",
	"grey85",
	"grey86",
	"grey87",
	"grey88",
	"grey89",
	"grey9",
	"grey90",
	"grey91",
	"grey92",
	"grey93",
	"grey94",
	"grey95",
	"grey96",
	"grey97",
	"grey98",
	"grey99",
	"honeydew",
	"honeydew1",
	"honeydew2",
	"honeydew3",
	"honeydew4",
	"hotPink",
	"hotPink1",
	"hotPink2",
	"hotPink3",
	"hotPink4",
	"indianRed",
	"indianRed1",
	"indianRed2",
	"indianRed3",
	"indianRed4",
	"ivory",
	"ivory1",
	"ivory2",
	"ivory3",
	"ivory4",
	"khaki",
	"khaki1",
	"khaki2",
	"khaki3",
	"khaki4",
	"lavender",
	"lavenderBlush",
	"lavenderBlush1",
	"lavenderBlush2",
	"lavenderBlush3",
	"lavenderBlush4",
	"lawnGreen",
	"lemonChiffon",
	"lemonChiffon1",
	"lemonChiffon2",
	"lemonChiffon3",
	"lemonChiffon4",
	"lightBlue",
	"lightBlue1",
	"lightBlue2",
	"lightBlue3",
	"lightBlue4",
	"lightCoral",
	"lightCyan",
	"lightCyan1",
	"lightCyan2",
	"lightCyan3",
	"lightCyan4",
	"lightGoldenrod",
	"lightGoldenrod1",
	"lightGoldenrod2",
	"lightGoldenrod3",
	"lightGoldenrod4",
	"lightGoldenrodYellow",
	"lightGray",
	"lightGreen",
	"lightGrey",
	"lightPink",
	"lightPink1",
	"lightPink2",
	"lightPink3",
	"lightPink4",
	"lightSalmon",
	"lightSalmon1",
	"lightSalmon2",
	"lightSalmon3",
	"lightSalmon4",
	"lightSeaGreen",
	"lightSkyBlue",
	"lightSkyBlue1",
	"lightSkyBlue2",
	"lightSkyBlue3",
	"lightSkyBlue4",
	"lightSlateBlue",
	"lightSlateGray",
	"lightSlateGrey",
	"lightSteelBlue",
	"lightSteelBlue1",
	"lightSteelBlue2",
	"lightSteelBlue3",
	"lightSteelBlue4",
	"lightYellow",
	"lightYellow1",
	"lightYellow2",
	"lightYellow3",
	"lightYellow4",
	"limeGreen",
	"linen",
	"magenta",
	"magenta1",
	"magenta2",
	"magenta3",
	"magenta4",
	"maroon",
	"maroon1",
	"maroon2",
	"maroon3",
	"maroon4",
	"mediumAquamarine",
	"mediumBlue",
	"mediumOrchid",
	"mediumOrchid1",
	"mediumOrchid2",
	"mediumOrchid3",
	"mediumOrchid4",
	"mediumPurple",
	"mediumPurple1",
	"mediumPurple2",
	"mediumPurple3",
	"mediumPurple4",
	"mediumSeaGreen",
	"mediumSlateBlue",
	"mediumSpringGreen",
	"mediumTurquoise",
	"mediumVioletRed",
	"midnightBlue",
	"mintCream",
	"mistyRose",
	"mistyRose1",
	"mistyRose2",
	"mistyRose3",
	"mistyRose4",
	"moccasin",
	"navajoWhite",
	"navajoWhite1",
	"navajoWhite2",
	"navajoWhite3",
	"navajoWhite4",
	"navy",
	"navyBlue",
	"oldLace",
	"oliveDrab",
	"oliveDrab1",
	"oliveDrab2",
	"oliveDrab3",
	"oliveDrab4",
	"orange",
	"orange1",
	"orange2",
	"orange3",
	"orange4",
	"orangeRed",
	"orangeRed1",
	"orangeRed2",
	"orangeRed3",
	"orangeRed4",
	"orchid",
	"orchid1",
	"orchid2",
	"orchid3",
	"orchid4",
	"paleGoldenrod",
	"paleGreen",
	"paleGreen1",
	"paleGreen2",
	"paleGreen3",
	"paleGreen4",
	"paleTurquoise",
	"paleTurquoise1",
	"paleTurquoise2",
	"paleTurquoise3",
	"paleTurquoise4",
	"paleVioletRed",
	"paleVioletRed1",
	"paleVioletRed2",
	"paleVioletRed3",
	"paleVioletRed4",
	"papayaWhip",
	"peachPuff",
	"peachPuff1",
	"peachPuff2",
	"peachPuff3",
	"peachPuff4",
	"peru",
	"pink",
	"pink1",
	"pink2",
	"pink3",
	"pink4",
	"plum",
	"plum1",
	"plum2",
	"plum3",
	"plum4",
	"powderBlue",
	"purple",
	"purple1",
	"purple2",
	"purple3",
	"purple4",
	"red",
	"red1",
	"red2",
	"red3",
	"red4",
	"rosyBrown",
	"rosyBrown1",
	"rosyBrown2",
	"rosyBrown3",
	"rosyBrown4",
	"royalBlue",
	"royalBlue1",
	"royalBlue2",
	"royalBlue3",
	"royalBlue4",
	"saddleBrown",
	"salmon",
	"salmon1",
	"salmon2",
	"salmon3",
	"salmon4",
	"sandyBrown",
	"seaGreen",
	"seaGreen1",
	"seaGreen2",
	"seaGreen3",
	"seaGreen4",
	"seashell",
	"seashell1",
	"seashell2",
	"seashell3",
	"seashell4",
	"sienna",
	"sienna1",
	"sienna2",
	"sienna3",
	"sienna4",
	"skyBlue",
	"skyBlue1",
	"skyBlue2",
	"skyBlue3",
	"skyBlue4",
	"slateBlue",
	"slateBlue1",
	"slateBlue2",
	"slateBlue3",
	"slateBlue4",
	"slateGray",
	"slateGray1",
	"slateGray2",
	"slateGray3",
	"slateGray4",
	"slateGrey",
	"snow",
	"snow1",
	"snow2",
	"snow3",
	"snow4",
	"springGreen",
	"springGreen1",
	"springGreen2",
	"springGreen3",
	"springGreen4",
	"steelBlue",
	"steelBlue1",
	"steelBlue2",
	"steelBlue3",
	"steelBlue4",
	"tan1",
	"tan2",
	"tan3",
	"tan4",
	"thistle",
	"thistle1",
	"thistle2",
	"thistle3",
	"thistle4",
	"tomato",
	"tomato1",
	"tomato2",
	"tomato3",
	"tomato4",
	"turquoise",
	"turquoise1",
	"turquoise2",
	"turquoise3",
	"turquoise4",
	"violet",
	"violetRed",
	"violetRed1",
	"violetRed2",
	"violetRed3",
	"violetRed4",
	"wheat",
	"wheat1",
	"wheat2",
	"wheat3",
	"wheat4",
	"white",
	"whiteSmoke",
	"yellow",
	"yellow1",
	"yellow2",
	"yellow3",
	"yellow4",
	"yellowGreen"
};

const static TColorGL rgbValue[RGBSIZE] =
{
	COLORaliceBlue,
	COLORantiqueWhite,
	COLORantiqueWhite1,
	COLORantiqueWhite2,
	COLORantiqueWhite3,
	COLORantiqueWhite4,
	COLORaquamarine,
	COLORaquamarine1,
	COLORaquamarine2,
	COLORaquamarine3,
	COLORaquamarine4,
	COLORazure,
	COLORazure1,
	COLORazure2,
	COLORazure3,
	COLORazure4,
	COLORbeige,
	COLORbisque,
	COLORbisque1,
	COLORbisque2,
	COLORbisque3,
	COLORbisque4,
	COLORblack,
	COLORblanchedAlmond,
	COLORblue,
	COLORblue1,
	COLORblue2,
	COLORblue3,
	COLORblue4,
	COLORblueViolet,
	COLORbrown,
	COLORbrown1,
	COLORbrown2,
	COLORbrown3,
	COLORbrown4,
	COLORburlywood,
	COLORburlywood1,
	COLORburlywood2,
	COLORburlywood3,
	COLORburlywood4,
	COLORcadetBlue,
	COLORcadetBlue1,
	COLORcadetBlue2,
	COLORcadetBlue3,
	COLORcadetBlue4,
	COLORchartreuse,
	COLORchartreuse1,
	COLORchartreuse2,
	COLORchartreuse3,
	COLORchartreuse4,
	COLORchocolate,
	COLORchocolate1,
	COLORchocolate2,
	COLORchocolate3,
	COLORchocolate4,
	COLORcoral,
	COLORcoral1,
	COLORcoral2,
	COLORcoral3,
	COLORcoral4,
	COLORcornflowerBlue,
	COLORcornsilk,
	COLORcornsilk1,
	COLORcornsilk2,
	COLORcornsilk3,
	COLORcornsilk4,
	COLORcyan,
	COLORcyan1,
	COLORcyan2,
	COLORcyan3,
	COLORcyan4,
	COLORdarkBlue,
	COLORdarkCyan,
	COLORdarkGoldenrod,
	COLORdarkGoldenrod1,
	COLORdarkGoldenrod2,
	COLORdarkGoldenrod3,
	COLORdarkGoldenrod4,
	COLORdarkGray,
	COLORdarkGreen,
	COLORdarkGrey,
	COLORdarkKhaki,
	COLORdarkMagenta,
	COLORdarkOliveGreen,
	COLORdarkOliveGreen1,
	COLORdarkOliveGreen2,
	COLORdarkOliveGreen3,
	COLORdarkOliveGreen4,
	COLORdarkOrange,
	COLORdarkOrange1,
	COLORdarkOrange2,
	COLORdarkOrange3,
	COLORdarkOrange4,
	COLORdarkOrchid,
	COLORdarkOrchid1,
	COLORdarkOrchid2,
	COLORdarkOrchid3,
	COLORdarkOrchid4,
	COLORdarkRed,
	COLORdarkSalmon,
	COLORdarkSeaGreen,
	COLORdarkSeaGreen1,
	COLORdarkSeaGreen2,
	COLORdarkSeaGreen3,
	COLORdarkSeaGreen4,
	COLORdarkSlateBlue,
	COLORdarkSlateGray,
	COLORdarkSlateGray1,
	COLORdarkSlateGray2,
	COLORdarkSlateGray3,
	COLORdarkSlateGray4,
	COLORdarkSlateGrey,
	COLORdarkTurquoise,
	COLORdarkViolet,
	COLORdeepPink,
	COLORdeepPink1,
	COLORdeepPink2,
	COLORdeepPink3,
	COLORdeepPink4,
	COLORdeepSkyBlue,
	COLORdeepSkyBlue1,
	COLORdeepSkyBlue2,
	COLORdeepSkyBlue3,
	COLORdeepSkyBlue4,
	COLORdimGray,
	COLORdimGrey,
	COLORdodgerBlue,
	COLORdodgerBlue1,
	COLORdodgerBlue2,
	COLORdodgerBlue3,
	COLORdodgerBlue4,
	COLORfirebrick,
	COLORfirebrick1,
	COLORfirebrick2,
	COLORfirebrick3,
	COLORfirebrick4,
	COLORfloralWhite,
	COLORforestGreen,
	COLORgainsboro,
	COLORghostWhite,
	COLORgold,
	COLORgold1,
	COLORgold2,
	COLORgold3,
	COLORgold4,
	COLORgoldenrod,
	COLORgoldenrod1,
	COLORgoldenrod2,
	COLORgoldenrod3,
	COLORgoldenrod4,
	COLORgray,
	COLORgray0,
	COLORgray1,
	COLORgray10,
	COLORgray100,
	COLORgray11,
	COLORgray12,
	COLORgray13,
	COLORgray14,
	COLORgray15,
	COLORgray16,
	COLORgray17,
	COLORgray18,
	COLORgray19,
	COLORgray2,
	COLORgray20,
	COLORgray21,
	COLORgray22,
	COLORgray23,
	COLORgray24,
	COLORgray25,
	COLORgray26,
	COLORgray27,
	COLORgray28,
	COLORgray29,
	COLORgray3,
	COLORgray30,
	COLORgray31,
	COLORgray32,
	COLORgray33,
	COLORgray34,
	COLORgray35,
	COLORgray36,
	COLORgray37,
	COLORgray38,
	COLORgray39,
	COLORgray4,
	COLORgray40,
	COLORgray41,
	COLORgray42,
	COLORgray43,
	COLORgray44,
	COLORgray45,
	COLORgray46,
	COLORgray47,
	COLORgray48,
	COLORgray49,
	COLORgray5,
	COLORgray50,
	COLORgray51,
	COLORgray52,
	COLORgray53,
	COLORgray54,
	COLORgray55,
	COLORgray56,
	COLORgray57,
	COLORgray58,
	COLORgray59,
	COLORgray6,
	COLORgray60,
	COLORgray61,
	COLORgray62,
	COLORgray63,
	COLORgray64,
	COLORgray65,
	COLORgray66,
	COLORgray67,
	COLORgray68,
	COLORgray69,
	COLORgray7,
	COLORgray70,
	COLORgray71,
	COLORgray72,
	COLORgray73,
	COLORgray74,
	COLORgray75,
	COLORgray76,
	COLORgray77,
	COLORgray78,
	COLORgray79,
	COLORgray8,
	COLORgray80,
	COLORgray81,
	COLORgray82,
	COLORgray83,
	COLORgray84,
	COLORgray85,
	COLORgray86,
	COLORgray87,
	COLORgray88,
	COLORgray89,
	COLORgray9,
	COLORgray90,
	COLORgray91,
	COLORgray92,
	COLORgray93,
	COLORgray94,
	COLORgray95,
	COLORgray96,
	COLORgray97,
	COLORgray98,
	COLORgray99,
	COLORgreen,
	COLORgreen1,
	COLORgreen2,
	COLORgreen3,
	COLORgreen4,
	COLORgreenYellow,
	COLORgrey,
	COLORgrey0,
	COLORgrey1,
	COLORgrey10,
	COLORgrey100,
	COLORgrey11,
	COLORgrey12,
	COLORgrey13,
	COLORgrey14,
	COLORgrey15,
	COLORgrey16,
	COLORgrey17,
	COLORgrey18,
	COLORgrey19,
	COLORgrey2,
	COLORgrey20,
	COLORgrey21,
	COLORgrey22,
	COLORgrey23,
	COLORgrey24,
	COLORgrey25,
	COLORgrey26,
	COLORgrey27,
	COLORgrey28,
	COLORgrey29,
	COLORgrey3,
	COLORgrey30,
	COLORgrey31,
	COLORgrey32,
	COLORgrey33,
	COLORgrey34,
	COLORgrey35,
	COLORgrey36,
	COLORgrey37,
	COLORgrey38,
	COLORgrey39,
	COLORgrey4,
	COLORgrey40,
	COLORgrey41,
	COLORgrey42,
	COLORgrey43,
	COLORgrey44,
	COLORgrey45,
	COLORgrey46,
	COLORgrey47,
	COLORgrey48,
	COLORgrey49,
	COLORgrey5,
	COLORgrey50,
	COLORgrey51,
	COLORgrey52,
	COLORgrey53,
	COLORgrey54,
	COLORgrey55,
	COLORgrey56,
	COLORgrey57,
	COLORgrey58,
	COLORgrey59,
	COLORgrey6,
	COLORgrey60,
	COLORgrey61,
	COLORgrey62,
	COLORgrey63,
	COLORgrey64,
	COLORgrey65,
	COLORgrey66,
	COLORgrey67,
	COLORgrey68,
	COLORgrey69,
	COLORgrey7,
	COLORgrey70,
	COLORgrey71,
	COLORgrey72,
	COLORgrey73,
	COLORgrey74,
	COLORgrey75,
	COLORgrey76,
	COLORgrey77,
	COLORgrey78,
	COLORgrey79,
	COLORgrey8,
	COLORgrey80,
	COLORgrey81,
	COLORgrey82,
	COLORgrey83,
	COLORgrey84,
	COLORgrey85,
	COLORgrey86,
	COLORgrey87,
	COLORgrey88,
	COLORgrey89,
	COLORgrey9,
	COLORgrey90,
	COLORgrey91,
	COLORgrey92,
	COLORgrey93,
	COLORgrey94,
	COLORgrey95,
	COLORgrey96,
	COLORgrey97,
	COLORgrey98,
	COLORgrey99,
	COLORhoneydew,
	COLORhoneydew1,
	COLORhoneydew2,
	COLORhoneydew3,
	COLORhoneydew4,
	COLORhotPink,
	COLORhotPink1,
	COLORhotPink2,
	COLORhotPink3,
	COLORhotPink4,
	COLORindianRed,
	COLORindianRed1,
	COLORindianRed2,
	COLORindianRed3,
	COLORindianRed4,
	COLORivory,
	COLORivory1,
	COLORivory2,
	COLORivory3,
	COLORivory4,
	COLORkhaki,
	COLORkhaki1,
	COLORkhaki2,
	COLORkhaki3,
	COLORkhaki4,
	COLORlavender,
	COLORlavenderBlush,
	COLORlavenderBlush1,
	COLORlavenderBlush2,
	COLORlavenderBlush3,
	COLORlavenderBlush4,
	COLORlawnGreen,
	COLORlemonChiffon,
	COLORlemonChiffon1,
	COLORlemonChiffon2,
	COLORlemonChiffon3,
	COLORlemonChiffon4,
	COLORlightBlue,
	COLORlightBlue1,
	COLORlightBlue2,
	COLORlightBlue3,
	COLORlightBlue4,
	COLORlightCoral,
	COLORlightCyan,
	COLORlightCyan1,
	COLORlightCyan2,
	COLORlightCyan3,
	COLORlightCyan4,
	COLORlightGoldenrod,
	COLORlightGoldenrod1,
	COLORlightGoldenrod2,
	COLORlightGoldenrod3,
	COLORlightGoldenrod4,
	COLORlightGoldenrodYellow,
	COLORlightGray,
	COLORlightGreen,
	COLORlightGrey,
	COLORlightPink,
	COLORlightPink1,
	COLORlightPink2,
	COLORlightPink3,
	COLORlightPink4,
	COLORlightSalmon,
	COLORlightSalmon1,
	COLORlightSalmon2,
	COLORlightSalmon3,
	COLORlightSalmon4,
	COLORlightSeaGreen,
	COLORlightSkyBlue,
	COLORlightSkyBlue1,
	COLORlightSkyBlue2,
	COLORlightSkyBlue3,
	COLORlightSkyBlue4,
	COLORlightSlateBlue,
	COLORlightSlateGray,
	COLORlightSlateGrey,
	COLORlightSteelBlue,
	COLORlightSteelBlue1,
	COLORlightSteelBlue2,
	COLORlightSteelBlue3,
	COLORlightSteelBlue4,
	COLORlightYellow,
	COLORlightYellow1,
	COLORlightYellow2,
	COLORlightYellow3,
	COLORlightYellow4,
	COLORlimeGreen,
	COLORlinen,
	COLORmagenta,
	COLORmagenta1,
	COLORmagenta2,
	COLORmagenta3,
	COLORmagenta4,
	COLORmaroon,
	COLORmaroon1,
	COLORmaroon2,
	COLORmaroon3,
	COLORmaroon4,
	COLORmediumAquamarine,
	COLORmediumBlue,
	COLORmediumOrchid,
	COLORmediumOrchid1,
	COLORmediumOrchid2,
	COLORmediumOrchid3,
	COLORmediumOrchid4,
	COLORmediumPurple,
	COLORmediumPurple1,
	COLORmediumPurple2,
	COLORmediumPurple3,
	COLORmediumPurple4,
	COLORmediumSeaGreen,
	COLORmediumSlateBlue,
	COLORmediumSpringGreen,
	COLORmediumTurquoise,
	COLORmediumVioletRed,
	COLORmidnightBlue,
	COLORmintCream,
	COLORmistyRose,
	COLORmistyRose1,
	COLORmistyRose2,
	COLORmistyRose3,
	COLORmistyRose4,
	COLORmoccasin,
	COLORnavajoWhite,
	COLORnavajoWhite1,
	COLORnavajoWhite2,
	COLORnavajoWhite3,
	COLORnavajoWhite4,
	COLORnavy,
	COLORnavyBlue,
	COLORoldLace,
	COLORoliveDrab,
	COLORoliveDrab1,
	COLORoliveDrab2,
	COLORoliveDrab3,
	COLORoliveDrab4,
	COLORorange,
	COLORorange1,
	COLORorange2,
	COLORorange3,
	COLORorange4,
	COLORorangeRed,
	COLORorangeRed1,
	COLORorangeRed2,
	COLORorangeRed3,
	COLORorangeRed4,
	COLORorchid,
	COLORorchid1,
	COLORorchid2,
	COLORorchid3,
	COLORorchid4,
	COLORpaleGoldenrod,
	COLORpaleGreen,
	COLORpaleGreen1,
	COLORpaleGreen2,
	COLORpaleGreen3,
	COLORpaleGreen4,
	COLORpaleTurquoise,
	COLORpaleTurquoise1,
	COLORpaleTurquoise2,
	COLORpaleTurquoise3,
	COLORpaleTurquoise4,
	COLORpaleVioletRed,
	COLORpaleVioletRed1,
	COLORpaleVioletRed2,
	COLORpaleVioletRed3,
	COLORpaleVioletRed4,
	COLORpapayaWhip,
	COLORpeachPuff,
	COLORpeachPuff1,
	COLORpeachPuff2,
	COLORpeachPuff3,
	COLORpeachPuff4,
	COLORperu,
	COLORpink,
	COLORpink1,
	COLORpink2,
	COLORpink3,
	COLORpink4,
	COLORplum,
	COLORplum1,
	COLORplum2,
	COLORplum3,
	COLORplum4,
	COLORpowderBlue,
	COLORpurple,
	COLORpurple1,
	COLORpurple2,
	COLORpurple3,
	COLORpurple4,
	COLORred,
	COLORred1,
	COLORred2,
	COLORred3,
	COLORred4,
	COLORrosyBrown,
	COLORrosyBrown1,
	COLORrosyBrown2,
	COLORrosyBrown3,
	COLORrosyBrown4,
	COLORroyalBlue,
	COLORroyalBlue1,
	COLORroyalBlue2,
	COLORroyalBlue3,
	COLORroyalBlue4,
	COLORsaddleBrown,
	COLORsalmon,
	COLORsalmon1,
	COLORsalmon2,
	COLORsalmon3,
	COLORsalmon4,
	COLORsandyBrown,
	COLORseaGreen,
	COLORseaGreen1,
	COLORseaGreen2,
	COLORseaGreen3,
	COLORseaGreen4,
	COLORseashell,
	COLORseashell1,
	COLORseashell2,
	COLORseashell3,
	COLORseashell4,
	COLORsienna,
	COLORsienna1,
	COLORsienna2,
	COLORsienna3,
	COLORsienna4,
	COLORskyBlue,
	COLORskyBlue1,
	COLORskyBlue2,
	COLORskyBlue3,
	COLORskyBlue4,
	COLORslateBlue,
	COLORslateBlue1,
	COLORslateBlue2,
	COLORslateBlue3,
	COLORslateBlue4,
	COLORslateGray,
	COLORslateGray1,
	COLORslateGray2,
	COLORslateGray3,
	COLORslateGray4,
	COLORslateGrey,
	COLORsnow,
	COLORsnow1,
	COLORsnow2,
	COLORsnow3,
	COLORsnow4,
	COLORspringGreen,
	COLORspringGreen1,
	COLORspringGreen2,
	COLORspringGreen3,
	COLORspringGreen4,
	COLORsteelBlue,
	COLORsteelBlue1,
	COLORsteelBlue2,
	COLORsteelBlue3,
	COLORsteelBlue4,
	COLORtan1,
	COLORtan2,
	COLORtan3,
	COLORtan4,
	COLORthistle,
	COLORthistle1,
	COLORthistle2,
	COLORthistle3,
	COLORthistle4,
	COLORtomato,
	COLORtomato1,
	COLORtomato2,
	COLORtomato3,
	COLORtomato4,
	COLORturquoise,
	COLORturquoise1,
	COLORturquoise2,
	COLORturquoise3,
	COLORturquoise4,
	COLORviolet,
	COLORvioletRed,
	COLORvioletRed1,
	COLORvioletRed2,
	COLORvioletRed3,
	COLORvioletRed4,
	COLORwheat,
	COLORwheat1,
	COLORwheat2,
	COLORwheat3,
	COLORwheat4,
	COLORwhite,
	COLORwhiteSmoke,
	COLORyellow,
	COLORyellow1,
	COLORyellow2,
	COLORyellow3,
	COLORyellow4,
	COLORyellowGreen
};



#endif