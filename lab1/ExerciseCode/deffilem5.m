%sets default properties for graphical objects
%(add deffile in the startup file matlabrc.m)

% TEXT PROPERTIES

	% Font name
%	set(0,'defaulttextfontname','Times')

	% Font size
	set(0,'defaulttextfontsize',16)

	% Color
%	set(0,'defaulttextcolor','k')

	% Font weight ( options:[ light | {normal} | demi | bold ] )
	% set(0,'defaulttextfontweight','normal')	

	% Font angle ( options:[ {normal} | italic | oblique ] )
%	set(0,'defaulttextfontangle','italic')	

	% Font strikethrough ( options:[ on | {off} ] )
%	set(0,'defaulttextfontstrikethrough','off')

	% Font underline ( options:[ on | {off} ] )
%	set(0,'defaulttextfontunderline','off')

	% Horizontal alignment ( options:[ {left} | center | right ] )
	set(0,'defaulttextHorizontalAlignment','left')

	% Vertical alignment ( options:[ top | cap | {middle} | baseline | bottom ] )
	set(0,'defaulttextverticalAlignment','middle')

	% Rotation
	set(0,'defaulttextrotation',0)

% LINE

	% Color
%	set(0,'defaultlinecolor','k')
	
	% Line style ( options:[ {-} | -- | : | -. | + | o | * | . | x ] )
	set(0,'defaultlinelinestyle','-')
	
	% Line width
	set(0,'defaultlinelinewidth',1.0)

	% Marker size
	set(0,'defaultlinemarkersize',10)




% FIGURE
	
	% Color
	set(0,'defaultfigurecolor','w')
	
	% Next plot ( options:[ new | {add} | replace ] )
	set(0,'defaultfigurenextplot','add')

	% Number title ( options:[ {on} | off ] )
	set(0,'defaultfigurenumbertitle','on')

	% Paper orientation ( options:[ {portrait} | landscape ] )
	set(0,'defaultfigurepaperorientation','portrait')

	% PaperType ( options:[ {usletter} | uslegal | a4letter ] )
	set(0,'defaultfigurepapertype','a4letter')

	% Pointer ( options:[ crosshair | {arrow} | watch | topl |
	%	    topr | botl | botr | circle | cross | fleur ] )
	set(0,'defaultfigurepointer','arrow')

	% Resize ( options:[ {on} | off ] )
	set(0,'defaultfigureresize','on')
	
	% Button down function
	set(0,'defaultfigurewindowbuttondownfcn','')
	
% AXES

	% Box ( options:[ on | {off} ] )
	set(0,'defaultaxesbox','off')
	
	% Color ( options:[ {none} ] -or- a ColorSpec. )
	set(0,'defaultaxescolor','w')

	% X Color 
%	set(0,'defaultaxesxcolor','k')
	
	% Y Color 
%	set(0,'defaultaxesycolor','k')

	% Font angle ( options:[ {normal} | italic | oblique ] )
	set(0,'defaultaxesfontangle','normal')

	% Font name
%	set(0,'defaultaxesfontname','times')

	% Font size (22 tidigare)
	set(0,'defaultaxesfontsize',16)

	% Font strikethrough ( options:[ on | {off} ] )
%	set(0,'defaultaxesfontstrikethrough','off')

	% Font underline ( options:[ on | {off} ] )
%	set(0,'defaultaxesfontunderline','off')

	% Font weight ( options:[ light | {normal} | demi | bold ] )
	set(0,'defaultaxesfontweight','normal')

	% Grid line style ( options:[ - | -- | {:} | -. ] )
	set(0,'defaultaxesgridlinestyle',':')

	% Line width
	set(0,'defaultaxeslinewidth',1.0)

	% Next plot ( options:[ new | add | {replace} ] )
	set(0,'defaultaxesnextplot','replace')

	% Tick length ( [2Dticklength  3Dticklength] )
	set(0,'defaultaxesticklength',[0.01 .025])

	% Tick direction ( options:[ {in} | out ]
	set(0,'defaultaxestickdir','in')


%	set(0,'DefaultAxesColorOrder',[0 0 0]);

	set(0,'DefaultAxesLineStyleOrder','-|--|:|-.');

