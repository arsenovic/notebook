# Alex Arsenovic's Notebook 

## Summary 
* This is a repo for my web-based notebook,   https://arsenovic.github.io/notebook/
* its a Work In Progress

## Why 
* I wanted to serve my jupyter notebooks online in an aesthetically pleasing way with minimal hassle. 

## Build Details 
For anyone who wants to re-use any part of this project, here are some details:
* this is a uv-based project
* The content  is stored in `content/` as jupytext-convertable markdown files. 
* The html is  built by `uv build.py`. Building details can be found in that file. 
* A github workflow allows the built html to be served on gh-pages. 



This repo is licensed under [CC BY-NC-ND](static/LICENSE.txt). 