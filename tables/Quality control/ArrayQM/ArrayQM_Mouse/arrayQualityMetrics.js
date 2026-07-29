// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, false, true, true, false, false, false, false, true, false, false, false, false, false ];
var arrayMetadata    = [ [ "1", "GSM3407171", "TriFLU-MF59_BL_1", "0", "Fluad_0_PBMC", "Mouse" ], [ "2", "GSM3407172", "TriFLU-MF59_BL_2", "0", "Fluad_0_PBMC", "Mouse" ], [ "3", "GSM3407173", "TriFLU-MF59_BL_3", "0", "Fluad_0_PBMC", "Mouse" ], [ "4", "GSM3407174", "TriFLU-MF59_BL_4", "0", "Fluad_0_PBMC", "Mouse" ], [ "5", "GSM3407175", "TriFLU-MF59_BL_5", "0", "Fluad_0_PBMC", "Mouse" ], [ "6", "GSM3407531", "TriFLU-MF59_BL_1", "1", "Fluad_24_PBMC", "Mouse" ], [ "7", "GSM3407532", "TriFLU-MF59_BL_2", "1", "Fluad_24_PBMC", "Mouse" ], [ "8", "GSM3407533", "TriFLU-MF59_BL_3", "1", "Fluad_24_PBMC", "Mouse" ], [ "9", "GSM3407534", "TriFLU-MF59_BL_4", "1", "Fluad_24_PBMC", "Mouse" ], [ "10", "GSM3407535", "TriFLU-MF59_BL_5", "1", "Fluad_24_PBMC", "Mouse" ], [ "11", "GSM3407770", "TriFLU-MF59_BL_1", "3", "Fluad_72_PBMC", "Mouse" ], [ "12", "GSM3407771", "TriFLU-MF59_BL_2", "3", "Fluad_72_PBMC", "Mouse" ], [ "13", "GSM3407772", "TriFLU-MF59_BL_3", "3", "Fluad_72_PBMC", "Mouse" ], [ "14", "GSM3407773", "TriFLU-MF59_BL_4", "3", "Fluad_72_PBMC", "Mouse" ], [ "15", "GSM3407774", "TriFLU-MF59_BL_5", "3", "Fluad_72_PBMC", "Mouse" ], [ "16", "GSM3407890", "TriFLU-MF59_BL_1", "7", "Fluad_168_PBMC", "Mouse" ], [ "17", "GSM3407891", "TriFLU-MF59_BL_2", "7", "Fluad_168_PBMC", "Mouse" ], [ "18", "GSM3407892", "TriFLU-MF59_BL_3", "7", "Fluad_168_PBMC", "Mouse" ], [ "19", "GSM3407893", "TriFLU-MF59_BL_4", "7", "Fluad_168_PBMC", "Mouse" ], [ "20", "GSM3407894", "TriFLU-MF59_BL_5", "7", "Fluad_168_PBMC", "Mouse" ] ];
var svgObjectNames   = [ "pca", "dens" ];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
    for(i=0; i<ssrules.length; i++) {
        if (ssrules[i].selectorText == (".aqm" + reportObjId)) {
		ssrules[i].style.cssText = cssText[0+status];
		break;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
