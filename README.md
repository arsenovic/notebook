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
* The html is  built by `uv run build.py`. Building details can be found in that file. 
* A github workflow allows the built html to be served on gh-pages. 


## TODO 
* make ganja.js animations work 
* use templating? could be used to support header/footer 
* automated way to datestamp. 
* katex vs mathjax?



This repo is licensed under [CC BY-NC-ND](static/LICENSE.txt). 