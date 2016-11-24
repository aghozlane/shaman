spincss <- "
  #plot-container {
z-index: 0;
position: relative;
}
#loading-spinner {
position: absolute;
left: 50%;
top: 50%;
z-index: -1;
margin-top: 33px;  /* half of the spinner's height */
margin-left: -33px; /* half of the spinner's width */
}
.recalculating {
z-index: -2;
background-color: #fff;
}
"

Errorcss <- "
.shiny-output-error { visibility: hidden; }
.shiny-output-error:before {
visibility: visible;
content: 'An error occurred. Please contact us at shaman@pasteur.fr'; }
}
"

appCSS <- "
#loading-content {
position: absolute;
background: #182b42;
opacity: 0.9;
z-index: 100;
left: 0;
right: 0;
top: 30px;
height: 100%;
text-align: center;
color: #FFFFFF;
}
#loading-content-bar {
position: absolute;
background: #182b42;
opacity: 0.9;
z-index: 100;
left: 0;
right: 0;
height: 100%;
text-align: center;
color: #FFFFFF;
}
"
