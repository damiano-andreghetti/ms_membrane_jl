<TeXmacs|2.1>

<style|generic>

<\body>
  <\itemize-minus>
    <item>le interazioni a lungo raggio mediate dalle fluttuazioni della
    membrana favoriscono la formazione di domini e aumentano l'efficienza del
    sorting

    <item>in subordine: questo potrebbe avvenire a causa di fenomeni di
    coalescenza di domini
  </itemize-minus>

  <paragraph|28/10/22>

  <\itemize-minus>
    <item><math|R=1/2>, vincolo è <math|d\<less\><sqrt|3>> <cite|KKN87|PB88>

    <item>mossa di Alexander, fare passaggio su tutti gli <math|N<rsub|e>>
    lati

    <item><with|font-shape|italic|accoppiare particelle e membrana>,
    curvatura sui vertici <with|color|blue|<math|\<longleftarrow\>>>

    <\itemize-minus>
      <item>meglio perché abbiamo reticolo approssimativamente esagonale
    </itemize-minus>

    <item>inserzione (dopo perché prima vogliamo capire se ci sono delle
    forze anche senza)

    <item>estrazione (dopo perché prima vogliamo capire se ci sono delle
    forze anche senza)

    <item><with|font-shape|italic|salvare configurazioni> per vedere durante
    la simulazione

    <item>scrivere Helfrich in vari modi discretizzati

    <item><math|k=0> <with|font-shape|italic|fa la fase tubulare>?

    <item>ottimizzazione Julia
  </itemize-minus>

  schema:

  <\itemize-minus>
    <item>fluttuazione membrana senza particelle, verificare dati di
    letteratura

    <\itemize-minus>
      <item>quali osservabili usano <cite|KKN87> per verificare la
      transizione di crumbling? la vediamo anche con la membrana chiusa?
      vedere <math|R<rsub|G><rsup|2>> vs <math|k> (rigidità) <cite|KKN87>

      <item>raggio quadratico di girazione medio (energia media) al variare
      di <math|k>
    </itemize-minus>
  </itemize-minus>

  Risultati:

  <\itemize-minus>
    <item><math|N<rsub|v>=642>, passo MC prende 4s
  </itemize-minus>

  Riferimento per superfici tethered:<nbsp><cite|KKN87>

  reply:

  https://journals-aps-org.unimib.idm.oclc.org/prl/abstract/10.1103/PhysRevLett.60.238

  absence of crumpling transition with self-avoidance:

  https://journals-aps-org.unimib.idm.oclc.org/pra/abstract/10.1103/PhysRevA.38.4943

  [qual'è la situazione? crumpling? roughening? verificare qualcosa che c'è
  in letteratura con nostro algoritmo?]

  spherical shells:

  https://journals-aps-org.unimib.idm.oclc.org/prx/abstract/10.1103/PhysRevX.7.011002

  thermal fluctuations increase the bending rigidity

  fluctuating shells under pressure

  https://www.pnas.org/doi/full/10.1073/pnas.1212268109

  review Gompper <cite|NPW04> verso pag. 350

  https://www-worldscientific-com.unimib.idm.oclc.org/doi/abs/10.1142/9789812565518_0012

  Membranes with Fluctuating Topology: Monte Carlo Simulations (Gompper,
  Kroll)

  membrane fluide che fanno transizione \Ptopologica\Q:

  https://iopscience-iop-org.unimib.idm.oclc.org/article/10.1088/0953-8984/9/42/001/pdf

  conformation of fluid membranes Gompper 1992:

  https://www-science-org.unimib.idm.oclc.org/doi/abs/10.1126/science.1546294

  <paragraph|17/10/22>

  <\itemize-minus>
    <item>perché non mettere la curvatura sui nodi? a quel punto basterebbe
    nella <math|H>, che è scritta come somma sui nodi <math|v>, inserire una
    <math|\<mu\>*\<delta\><rsub|v>>

    <item>verificare se possibile che la curvatura calcolata con lo shape
    operator tiene conto correttamente anche dei contributi di curvatura dai
    lati

    <item>che contributo dà un termine <math|H<rsub|0>=2/R<rsub|0>> costante
    dentro a <math|<around*|(|H<rsub|<inactive*|c>>-H<rsub|0>|)><rsup|2>>

    <item>esiste qualche modo di evitare il controllo degli overlap?

    <item>la <math|\<ell\>> nella tesi di Riccardo Rossetti (lunghezza del
    filo) è contata da superficie a superficie?

    <item>implementare: link flip, inserzione, estrazione

    <item>vedere se due particelle (o due domini) si attraggono

    <item>vedere in letteratura analoghe prove dell'interazione a lungo
    raggio tra le inclusioni

    <item>se si attraggono producono una specie di <math|g<rsub|eff>\<gtr\>g>

    <item>se si attraggono potrebbe vedersi una distribuzione non uniforme di
    domini

    <item>quali sono le caratteristiche della distribuzione di domini e di
    distanze tra i domini?

    <item>confrontare distribuzione distanze tra i domini con membrana che
    fluttua/non fluttua

    <item>scrivere riassunto

    <item>guardare crumbling

    <item>correzione per diffusione da triangoli non equilateri (rate
    proporzionale al lato?)
  </itemize-minus>

  <\eqnarray*>
    <tformat|<table|<row|<cell|k,\<mu\>,\<phi\>,g>|<cell|>|<cell|>>>>
  </eqnarray*>

  Possibili esiti desiderabili:

  <\itemize-minus>
    <item>evidenza di attrazione tra le inclusioni (di base, vedere cosa c'è
    in letteratura)\ 

    <item>distillazione è favorita/sfavorita da interazione mediata da
    membrana

    <item>distribuzione distanze dei domini non banale
  </itemize-minus>

  <\bibliography|bib|alpha|/home/andrea/archivio/Archivio>
    <\bib-list|KKN87>
      <bibitem*|KKN87><label|bib-KKN87>Yacov Kantor, Mehran Kardar, and
      David<nbsp>R Nelson. <newblock>Tethered surfaces: Statics and dynamics.
      <newblock><with|font-shape|italic|Physical Review A>, 35(7):3056, 1987.

      <bibitem*|PB88><label|bib-PB88>Michael Plischke and David Boal.
      <newblock>Absence of a crumpling transition in strongly self-avoiding
      tethered membranes. <newblock><with|font-shape|italic|<pra>>,
      38(9):4943\U4945, November 1988.
    </bib-list>
  </bibliography>
</body>

<\initial>
  <\collection>
    <associate|page-medium|paper>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-2|<tuple|2|2>>
    <associate|auto-3|<tuple|<with|mode|<quote|math>|<rigid|->>|2>>
    <associate|bib-KKN87|<tuple|KKN87|2>>
    <associate|bib-PB88|<tuple|PB88|2>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      KKN87

      PB88

      KKN87

      KKN87

      KKN87

      NPW04
    </associate>
    <\associate|toc>
      <with|par-left|<quote|4tab>|28/10/22
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|17/10/22
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.15fn>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>