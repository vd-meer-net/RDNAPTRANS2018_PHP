# rdnaptrans2018_php

# Inleiding

<p>Deze zeer eenvoudige PHP-pagina <b>rdnaptrans.php</b> is een PHP-implementatie die de geografische Nederlandse <b>RD NAP</b> (<b>R</b>ijks<b>D</b>riehoeksmeting en <b>N</b>ormaal
<b>A</b>msterdam <b>P</b>eil) coördinaten omzet, transformeert, naar <b>ETRS89</b> (<b>E</b>uropean <b>T</b>errestrial <b>R</b>eference
<b>S</b>ystem 1989), of andersom. De code is gecertificeerd en mag het handelsmerk <b>RDNAPTRANS™2018</b> voeren. Dit betekent dat deze transformaties correct zijn als er juist gebruik wordt gemaakt van <b>rdnaptrans2018.php</b>. <b>RD NAP</b> eenheid is in meters, <b>ETRS89</b> is in graden en meters (de hoogte).</p>

<p><b>RDNAPTRANS™2018</b> compliant code transformeert elk punt (binnen <strong>RD NAP</strong> en <strong>ETRS89</strong> domein), welke plek op aarde dan ook. Maar buiten de zogenaamde grids kan de afwijking groot zijn en klopt er niets meer van. Dat is volkomen correct gedrag. Sommige implementaties geven dan een waarschuwing dat je transformeert met waarden die buiten het grid liggen. Deze code geeft die waarschuwing (nog) niet.</p>

<p><strong><u>Let op:</u></strong> De twee validatiebestanden nodig voor de zelfservice certificering zijn aan verandering onderhevig. De punten die er in staan veranderen. Beter is om ze van de <strong>NSGI</strong>-website zelf te halen. Dan heb je altijd de laatste versie. Voor de werking van PHP-pagina heb je ze niet nodig.</p>

<p>Tot zover de inleiding. Het is ook allemaal te vinden op <a href="http://www.nsgi.nl">www.nsgi.nl</a>.</p>

# De Bestanden
<p>De PHP-code is opgebouwd door elke stap na te rekenen met de <a href="https://github.com/FVellinga/gm_rdnaptrans2018/blob/main/README.md">SAS-implementatie van RDNAPTRANS™2018</a>.</p>

<p>De index.php toont 5 keuzes:</p>
<ul>
<li>Convert 1 demo-coordinate to ETRS89</li>
<li>Convert 1 demo-coordinate to RD</li>
<li>Perform selfvalidation</li>
<li>Create validationfile for ETRS2RD</li>
<li>Create validationfile for RD2ETRS</li>
</ul>

<p>De eerste twee zetten een coördinaat in de index.php om van RD naar ETRS89 of v.v.</p>
<p>De derde keuze laadt het bestand in met zelfvalidatie-punten en zet de coördinaten in beide richtingen. Indien het antwoord afwijkt van de invoer, wordt deze rood gekleurd. Dit is geen probleem, omdat de afwijking binnen de marges valt.</p>
<p>De vierde en vijfde keuze lezen de valiatiebestanden in en converteren de coördinaten. Het resultaat wordt aangeboden als download.</p>

<p>Bijgevoegd zijn er twee gridfiles (tekstbestanden). Deze worden ingeladen en gebruikt tijens de conversie. Dit betreffen rdcorr2018.txt en nlgeo2018.txt. Deze bestanden zijn te downloaden op de eereder genoemde websites van NGIS. 
Optioneel zijn het zogenaamde zelfvalidatiebestand (Z001_ETRS89andRDNAP.txt) en de twee certificatie validatiebestanden (002_ETRS90.txt en 002_RDNAP.txt). Handig, want hiermee toon je aan dat deze code werkt.</p>

# Het gebruik

<p>De index.php is niet nodig voor de berekening, die vormt alleen maar een eenvoudige schil.</p>
<p>Alle berekeningen worden uitgevoerd in rdnaptrans.php. Dit PHP-bestand bevat alle functies die nodig zijn voor de berekening. De enige twee functies die je hoeft aan te roepen zijn:</p>
<ul>
<li>etrs2rd()</li>
<li>rd2etrs2()</li>
</ul>

# etrs2rd()
<p>Deze functie vraagt maar 1 parameter en dat is een array met daarin de latitude en longitude )en eventueel hoogte) van het prunt dat omgerekend moet worden naar RD.</p>
<code>$rd=etrs2rd(array('lat'=>53.473095072,'lon'=>6.886681267,'h'=>290.4741));</code>

<p>Na de berekening bevat $rd een array, met daarin een x, een y en een h.</p>

# rd2etrs()
<p>Deze functie vraagt maar 1 parameter en dat is een array met daarin de x en de y ) en eventueel hoogte) van het prunt dat omgerekend moet worden naar ETRS.</p>
<code>$etrs=rd2etrs(array('x'=>108360.8790,'y'=>415757.2745,'h'=>258.0057))</code>

<p>Na de berekening bevat $etrs een array, met daarin een lat, een lon en een h.</p>

