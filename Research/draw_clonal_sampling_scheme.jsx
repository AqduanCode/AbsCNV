#target illustrator
(function () {
  var W = 496.06; // 175 mm at 72 pt/in
  var H = 360;
  var outFolder=Folder.selectDialog('Choose a new output folder for the clonal sampling scheme');
  if(!outFolder) return;
  var outAI=new File(outFolder.fsName+'/clonal_sampling_scheme.ai');
  var outPNG=new File(outFolder.fsName+'/clonal_sampling_scheme');
  var outLog=new File(outFolder.fsName+'/clonal_sampling_scheme_export.txt');
  function color(hex) {
    var c = new RGBColor();
    c.red = parseInt(hex.substr(1, 2), 16);
    c.green = parseInt(hex.substr(3, 2), 16);
    c.blue = parseInt(hex.substr(5, 2), 16);
    return c;
  }
  var ink = color('#263247'), muted = color('#596779'), edge = color('#C7D0DA');
  var purple = color('#A24E83'), teal = color('#287F88'), blue = color('#4268A7');
  var paleP = color('#F5EBF1'), paleT = color('#EAF4F4'), paleB = color('#EDF2F9');
  var white = color('#FFFFFF');
  var doc = app.documents.add(DocumentColorSpace.RGB, W, H);
  doc.layers[0].name = 'Connectors';
  var links = doc.layers[0];
  var cells = doc.layers.add(); cells.name = 'Editable cells and samples';
  var labels = doc.layers.add(); labels.name = 'Editable text';
  var regular = app.textFonts.getByName('ArialMT');
  var bold = app.textFonts.getByName('Arial-BoldMT');
  function path(layer, points, stroke, sw, fill, name) {
    var p = layer.pathItems.add();
    var arr = [];
    for (var i = 0; i < points.length; i++) arr.push([points[i][0], H - points[i][1]]);
    p.setEntirePath(arr);
    p.closed = !!fill;
    p.stroked = !!stroke;
    if (stroke) { p.strokeColor = stroke; p.strokeWidth = sw || 0.8; }
    p.filled = !!fill;
    if (fill) p.fillColor = fill;
    p.name = name || 'vector path';
    return p;
  }
  function line(layer, points, stroke, sw, name) { return path(layer, points, stroke, sw, null, name); }
  function arrow(layer, x1, y1, x2, y2, stroke, sw, name) {
    line(layer, [[x1,y1],[x2,y2]], stroke, sw, name || 'arrow shaft');
    var dx=x2-x1, dy=y2-y1, len=Math.sqrt(dx*dx+dy*dy), ux=dx/len, uy=dy/len;
    var px=-uy, py=ux, s=3.4;
    path(layer,[[x2,y2],[x2-ux*s+px*s*.55,y2-uy*s+py*s*.55],[x2-ux*s-px*s*.55,y2-uy*s-py*s*.55]],null,0,stroke,'arrowhead');
  }
  function circle(layer, cx, cy, r, fill, stroke, sw, name) {
    var p = layer.pathItems.ellipse(H-(cy-r),cx-r,2*r,2*r);
    p.filled=!!fill; if(fill)p.fillColor=fill;
    p.stroked=!!stroke; if(stroke){p.strokeColor=stroke;p.strokeWidth=sw||.8;}
    p.name=name||'editable cell circle';return p;
  }
  function roundRect(layer,x,y,w,h,rad,fill,stroke,sw,name){
    var p=layer.pathItems.roundedRectangle(H-y,x,w,h,rad,rad);
    p.filled=!!fill;if(fill)p.fillColor=fill;
    p.stroked=!!stroke;if(stroke){p.strokeColor=stroke;p.strokeWidth=sw||.7;}
    p.name=name||'editable rounded rectangle';return p;
  }
  function text(s,x,y,size,fill,isBold,center,name) {
    var t=labels.textFrames.add();t.contents=s;t.position=[x,H-y];
    var a=t.textRange.characterAttributes;
    a.textFont=isBold?bold:regular;a.size=size;a.fillColor=fill;
    if(center)t.left=x-t.width/2;
    t.name=name||s;return t;
  }
  function singleCell(x,y,tone,r,name) {
    circle(cells,x,y,r,white,tone,1.0,name+' outer');
    circle(cells,x,y,r*.49,tone,null,0,name+' nucleus');
  }
  function culture(x,y,tone,r,name) {
    var pts=[[-r*.62,-r*.22],[0,-r*.76],[r*.62,-r*.22],[-r*.43,r*.59],[r*.42,r*.58]];
    for(var i=0;i<pts.length;i++)singleCell(x+pts[i][0],y+pts[i][1],tone,r*.42,name+' cell '+(i+1));
  }
  function sample(x,y,w,tone,light,name) {
    roundRect(cells,x-w/2,y-6,w,12,4,light,tone,.55,name+' sample box');
    text('DNA  |  RNA',x,y-9,6.2,tone,true,true,name+' sample label');
  }
  // Title and methodological scale.
  text('U251 clonal lineage and sampling',18,23,11.4,ink,true,false,'figure title');
  text('U1 parent  /  3 P1 clones  /  2 P2 descendants per P1 clone  /  10 samples',18,38,7.25,muted,false,false,'sample count');

  // Parent/U1: initiating cell, 30-day expansion, DNA/RNA sampling.
  singleCell(248,61,purple,5.7,'parent U1');
  text('U1 parental culture',248,70,8.2,ink,true,true,'U1 label');
  arrow(links,248,90,248,101,purple,1.2,'first 30-day expansion');
  text('~30 d',267,95,6.8,muted,false,false,'first duration');
  culture(248,113,purple,9,'expanded parent');
  line(links,[[257,113],[300,113]],purple,.85,'parent sample branch');
  sample(336,113,68,purple,paleP,'U1 DNA RNA');

  // First single-cell isolation into three P1 clones.
  text('Single-cell isolation',270,123,7.2,purple,true,false,'first cloning label');
  line(links,[[248,126],[248,141],[87,141],[409,141]],teal,1.0,'U1 to P1 lineage');
  var p1=[87,248,409];
  for(var i=0;i<3;i++){
    var x=p1[i];
    line(links,[[x,141],[x,152]],teal,1.0,'P1 '+(i+1)+' branch');
    singleCell(x,160,teal,4.6,'P1 '+(i+1)+' seed');
    arrow(links,x,166,x,194,teal,1.0,'P1 '+(i+1)+' expansion');
    text('~30 d',x+9,183,6.3,muted,false,false,'P1 '+(i+1)+' duration');
    culture(x,205,teal,8,'P1 '+(i+1)+' culture');
    text('P1 clone '+(i+1),x+12,218,7.35,ink,true,false,'P1 '+(i+1)+' label');
    line(links,[[x+8,205],[x+27,205]],teal,.75,'P1 '+(i+1)+' sample branch');
    sample(x+52,205,50,teal,paleT,'P1 '+(i+1)+' DNA RNA');
  }
  // Second single-cell isolation, two descendants from each P1 culture.
  text('Single-cell isolation',264,233,6.6,blue,true,false,'second cloning label');
  var p2=[];
  for(var j=0;j<3;j++){
    var px=p1[j], left=px-37, right=px+37;
    line(links,[[px,213],[px,253],[left,253],[right,253]],blue,.95,'P1 '+(j+1)+' to P2 branch');
    p2.push(left);p2.push(right);
  }
  for(var k=0;k<6;k++){
    var xx=p2[k];
    line(links,[[xx,253],[xx,265]],blue,.9,'P2 '+(k+1)+' branch');
    singleCell(xx,270,blue,3.8,'P2 '+(k+1)+' seed');
    arrow(links,xx,275,xx,295,blue,.9,'P2 '+(k+1)+' expansion');
    text('~30 d',xx+5,288,5.6,muted,false,false,'P2 '+(k+1)+' duration');
    culture(xx,307,blue,6.7,'P2 '+(k+1)+' culture');
    text('P2 '+(k+1),xx+11,304,6.8,ink,true,false,'P2 '+(k+1)+' label');
    sample(xx,335,54,blue,paleB,'P2 '+(k+1)+' DNA RNA');
  }

  // Save the native document and an Illustrator-rendered visual check.
  var opts=new IllustratorSaveOptions();opts.pdfCompatible=true;opts.compressed=true;
  doc.saveAs(outAI,opts);
  var exp=new ExportOptionsPNG24();exp.antiAliasing=true;exp.artBoardClipping=true;
  exp.transparency=false;exp.horizontalScale=250;exp.verticalScale=250;
  doc.exportFile(outPNG,ExportType.PNG24,exp);
  outLog.encoding='UTF-8';outLog.open('w');
  outLog.writeln('Saved AI: '+outAI.fsName);
  outLog.writeln('Layers: '+doc.layers.length);
  outLog.writeln('Text frames: '+doc.textFrames.length);
  outLog.writeln('Path items: '+doc.pathItems.length);
  outLog.writeln('Placed images: '+doc.placedItems.length);
  outLog.writeln('Raster images: '+doc.rasterItems.length);
  outLog.writeln('Artboard pt: '+W+' x '+H);
  outLog.close();
  app.redraw();
})();
