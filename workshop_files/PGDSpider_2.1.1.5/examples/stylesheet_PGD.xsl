<?xml version="1.0"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<!-- stylesheet of the data format version 6 -->

 <xsl:template match="/">
  <html>
   <head>
    <title/>

<!-- xxxx  style definitions  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx -->

    <style type="text/css"> 
      h1 {font-family:Arial; font-size:18px; color:#00c;}
      h3 {font-family:Arial; font-size:16px; color:#00c}
      p {font-family:Arial; font-size:12px; font-weight:bold;}
      span {font-family:Arial; font-size:12px; font-weight:normal;}
      table {font-family:Arial; font-size:12px; border=0; border-collapse:collapse;}
      th {text-align:left; border:1px solid gray; padding-left:5px; padding-right:5px}
      th.without {text-align:left; border:0px; padding-left:5px; padding-right:5px}
      th.withoutInd {color:#00c; text-align:left; border:0px; padding-left:5px; padding-right:5px}
      th.matrix {text-align:left; vertical-align:text-top; border:1px solid gray; padding-left:0px; padding-right:5px}
      td {text-align:left; border:1px solid gray; padding-left:5px; padding-right:5px}
      td.without {text-align:left; border:0px; padding-left:5px; padding-right:5px}
      td.matrix {text-align:left; border:1px solid gray; padding-left:5px; padding-right:5px;}
      td.genotype {font-family:Lucida Console; font-size:12px; text-align:left; border:1px solid gray; padding-left:5px; padding-right:5px}
      td.genotype2 {font-family:Lucida Console; font-size:12px; text-align:left; border:0px; padding-left:5px; padding-right:5px}    
    </style>

   </head>

   <body>
    <xsl:apply-templates/>
   </body>
  </html>
 </xsl:template>



<!-- xxxx  header  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx -->
 
 <xsl:template match="header">
    <h1>
      <xsl:value-of select="@title"/>
    </h1>

    <table>
      <xsl:if test="organism">
        <tr>
          <th class="without"> Organism:</th> 
          <td class="without">
            <xsl:value-of select="organism"/>
          </td>
        </tr>
      </xsl:if>
      
      <xsl:if test="numPop">
        <tr>
          <th class="without">Number of populations:</th> 
          <td class="without">
            <xsl:value-of select="numPop"/>
          </td>
        </tr>
      </xsl:if>

      <xsl:if test="ploidy">
        <tr>
          <th class="without"> Ploidy: </th>
          <td class="without">
            <xsl:value-of select="ploidy"/> 
          </td>
        </tr>
      </xsl:if>

      <xsl:if test="missing">
        <tr>
          <th class="without"> Symbol for missing data: </th>
          <td class="without">
            <xsl:value-of select="missing"/>
          </td>
        </tr>
      </xsl:if>

      <xsl:if test="gap">
        <tr>
          <th class="without"> Symbol for gaps: </th>
          <td class="without">
            <xsl:value-of select="gap"/>
          </td> 
        </tr> 
      </xsl:if>      

      <xsl:if test="gameticPhase">
        <tr>
          <th class="without"> Gametic phase: </th>
          <td class="without"> known </td>
        </tr>
      </xsl:if>

      <xsl:if test="recessiveData">
        <tr>
          <th class="without"> recessive data: </th>
          <td class="without"> yes </td>
        </tr>
      </xsl:if> 
     
   </table>
 </xsl:template>







<!-- xxxx  dataDescription  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx -->

<xsl:template match="dataDescription">
  <h3> 
   <br />
   Data description: 
  </h3>
  <p> 
    <xsl:if test="numLoci">
      Loci number: 
      <span> 
        <xsl:value-of select="numLoci" />
        <br />
      </span>
   </xsl:if>

   <xsl:if test="dataType">
     Data type:
     <span> 
       <xsl:value-of select="dataType" />
       <br />
     </span>
   </xsl:if>
  </p>
  
  <table>
   <xsl:for-each select="locus">
    <tr>    
      <th class="withoutInd">
        id: <xsl:value-of select="@id" />
      </th>
    </tr>

    <xsl:if test="locusDataType">
      <tr>
        <td class="without"> data type: </td>
        <td class="without">
          <xsl:value-of select="locusDataType" />
        </td>
      </tr>
    </xsl:if>

    <xsl:if test="locusChromosome">
      <tr>
        <td class="without"> on chromosome: </td>
        <td class="without">
          <xsl:value-of select="locusChromosome" />
        </td>
      </tr>
    </xsl:if>

    <xsl:if test="locusLocation">
      <tr>
        <td class="without"> location: </td>
        <td class="without">
          <xsl:value-of select="locusLocation" />
        </td>
      </tr>
    </xsl:if>

    <xsl:if test="locusGenic">
      <tr>
        <td class="without"> locus is genic: </td>
        <td class="without">
          <xsl:value-of select="locusGenic" />
        </td>
      </tr>
    </xsl:if>

    <xsl:if test="locusLength">
      <tr>        
        <td class="without"> length: </td>
        <td class="without">
          <xsl:value-of select="locusLength" />
        </td>     
      </tr>
    </xsl:if>
    
    <xsl:if test="locusAncestralState">
      <tr>
        <td class="without"> ancestral state: </td>
        <td class="without">
          <xsl:value-of select="locusAncestralState" />
        </td>
      </tr>
    </xsl:if>

    <xsl:if test="locusLinks">
      <tr>         
        <td class="without"> links: </td>
        <td class="without">
          <xsl:value-of select="locusLinks" />
        </td>     
      </tr>
    </xsl:if>

    <xsl:if test="locusComments">
      <tr>
        <td class="without"> comments: </td>
        <td class="without">
          <xsl:value-of select="locusComments" />
        </td>
      </tr>
    </xsl:if>

    <tr>  
      <td class="without">
        <br /> 
      </td>
    </tr>

   </xsl:for-each>
  </table>
</xsl:template>





  

<!-- xxxx  population  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx -->

 <xsl:template match="population">


<!-- xxxx loci begin for alignment (variable setting) xxxxxxxxxx -->

     <xsl:variable name="min">
       <xsl:for-each select="ind/read/start">
         <xsl:sort select="." />
         <xsl:if test="position() = 1">
           <xsl:value-of select="."/>
         </xsl:if> 
       </xsl:for-each>
     </xsl:variable>

<!-- xxxxxxxxxxxxxx -->


  <h3> 
   <br />
   Population: <xsl:value-of select="@name"/>
  </h3>
  <p> Population size:
   <span> 
    <xsl:value-of select="popSize" />
    <br />
   </span>  

   <xsl:if test="popGeogCoord">
     Geographic coordination:
     <span>
       <xsl:value-of select="popGeogCoord"/>
     </span>
     <br /> 
   </xsl:if> 

   <xsl:if test="popLingGroup">
     Linguistic group:
     <span>
       <xsl:value-of select="popLingGroup"/>
     </span>
     <br />
   </xsl:if>


   <xsl:if test="popPloidy">
     Ploidy:
     <span>
       <xsl:value-of select="popPloidy"/>
     </span>
     <br />
   </xsl:if>


   <xsl:if test="popLoci">
     Loci:
     <span>
       <xsl:value-of select="popLoci" />
     </span>
     <br />
     <xsl:if test="ind/read/start">
       Locus start:
       <span>
         <xsl:value-of select="$min" />
       </span>
       <br />
     </xsl:if>
   </xsl:if>  
  </p>


<!-- xxxx variable if indFreq tag exist or not xxxxxxxx -->

   <xsl:variable name="Freq">
     <xsl:choose>
       <xsl:when test="ind/indFreq">
         1
       </xsl:when>
       <xsl:otherwise>
         0
       </xsl:otherwise>
     </xsl:choose>
   </xsl:variable>

<!-- xxxxxxxxxxxxxx -->


   <table>
     <xsl:apply-templates select="ind">
       <xsl:with-param name="minimum" select="$min + 1" /> 
       <xsl:with-param name="frequency" select="$Freq" />
     </xsl:apply-templates>
   </table>
 </xsl:template> 






<!-- xxxx  ind  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx -->

 <xsl:template match="ind">
  <xsl:param name="minimum" /> 
  <xsl:param name="frequency" />



<!-- xxx align over one individual (different loci) xxxxxxxxxxxxxxxxxxxx-->

   <xsl:choose>
     <xsl:when test="indLoci">
       <tr>
         <th class="withoutInd"> Individual <xsl:value-of select="@name" />: </th>
       </tr>

       <xsl:if test="indGeogCoord">
         <tr>
           <td class="without"> geographic coordination: </td>
           <td class="without">
             <xsl:value-of select="indGeogCoord" />
           </td>
         </tr>
       </xsl:if>

       <xsl:if test="indLingGroup">
         <tr>
           <td class="without"> linguistic group: </td>
           <td class="without"> 
             <xsl:value-of select="indLingGroup" />
           </td>
         </tr>
       </xsl:if>
       
       <xsl:if test="indLoci">
         <tr>
           <td class="without"> Loci: </td>
           <td class="without"> 
             <xsl:value-of select="indLoci" />
           </td>
         </tr>
       </xsl:if>



<!-- xxxx loci start for alignment (variable setting) xxxxxxx -->

       <xsl:variable name="min">
         <xsl:for-each select="read/start">
           <xsl:sort select="." />
           <xsl:if test="position() = 1">
             <xsl:value-of select="."/>
           </xsl:if> 
         </xsl:for-each>
       </xsl:variable>

<!-- xxxxxxxxxxxxxxxxxxxx -->


       <xsl:if test="$min &gt; 0">
         <tr>
           <td class="without"> Locus start: </td>
           <td class="without">
             <xsl:value-of select="$min" />
           </td>
         </tr>
       </xsl:if>

       <xsl:if test="indFreq">
         <tr>
           <td class="without"> Frequency of this individual: </td>
           <td class="without">
             <xsl:value-of select="indFreq" />
           </td>
         </tr>
       </xsl:if>

       <xsl:if test="indPloidy">
         <tr>
           <td class="without"> Ploidy: </td>
           <td class="without">
             <xsl:value-of select="indPloidy" />
           </td>
         </tr>
       </xsl:if>


       <tr>
         <td class="without"> Data: </td>
         


<!-- xxxx data alignment xxxxxxxxxxxxxxxxx -->
         <td class="genotype2">
           <xsl:choose>
             <xsl:when test="read">
               <xsl:apply-templates select="read">
                 <xsl:with-param name="minimum" select="$min + 1" />
               </xsl:apply-templates>
             </xsl:when>

             <!-- aligned data -->
             <xsl:otherwise>
               <xsl:for-each select="data">
                 <xsl:value-of select="." /> 
                 <br />
               </xsl:for-each>
             </xsl:otherwise>
           </xsl:choose>
         </td>
       </tr>

       <tr>  
         <td class="without">
           <br /> 
         </td>
       </tr>         
     
     </xsl:when>




<!-- xxx align over one population (ind same loci) xxxxxxxxxxxxxxxxxxxxxx-->

     <xsl:otherwise>

       <xsl:if test="position()=1">
         <tr>
           <th>ind</th>

           <xsl:if test="indGeogCoord">
             <th> geog. coord. </th>
           </xsl:if>

           <xsl:if test="indLingGroup">
             <th> ling. group </th>
           </xsl:if>

           <xsl:if test="$frequency = 1">
             <th> freq </th>
           </xsl:if>

           <xsl:if test="indPloidy">
             <th> ploidy </th>
           </xsl:if>
      
           <xsl:if test="data">
             <th>
               data 

               <!-- xxx test if the popLoci siblings before the current node exist xxx --> 
               <xsl:if test="preceding-sibling::popLoci">
                 <span>
                   <br />
                   (<xsl:value-of select="preceding-sibling::popLoci" />)
                 </span> 
               </xsl:if>
             </th> 
           </xsl:if>       

           <xsl:if test="read/data">
             <th> data </th>
           </xsl:if>

         </tr>
       </xsl:if>

       <tr>  
         <td> <xsl:value-of select="@name"/> </td>

         <xsl:if test="indGeogCoord">
           <td> 
             <xsl:value-of select="indGeogCoord" />
           </td>     
         </xsl:if>

         <xsl:if test="indLingGroup">
           <td> 
             <xsl:value-of select="indLingGroup" />
          </td>   
         </xsl:if>

         <xsl:if test="$frequency = 1">
           <td>
             <xsl:choose>
               <xsl:when test="indFreq"> 
                 <xsl:value-of select="indFreq" />
               </xsl:when>
               <xsl:otherwise>
                 1
               </xsl:otherwise>
             </xsl:choose>
          </td>   
         </xsl:if>

         <xsl:if test="indPloidy">
           <td> 
             <xsl:value-of select="indPloidy" />
           </td>   
         </xsl:if>
         
         <xsl:if test="data">
           <td class="genotype">    
             <xsl:for-each select="data">
               <xsl:value-of select="." /> 
               <br />
             </xsl:for-each>
           </td>
         </xsl:if>

         <xsl:if test="read">
           <td class="genotype">
             <xsl:apply-templates select="read">
               <xsl:with-param name="minimum" select="$minimum" />
             </xsl:apply-templates>
           </td>
         </xsl:if>

       </tr>

     </xsl:otherwise>
   </xsl:choose>
 </xsl:template>



<!-- xxxx template read xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx -->

 <xsl:template match="read"> 
  <xsl:param name="minimum" />
 
       <xsl:call-template name="align">
         <xsl:with-param name="min" select="$minimum" />
         <xsl:with-param name="beginning" select="start" />
       </xsl:call-template>

       <xsl:value-of select="data" />
       <br />

 </xsl:template>



<!-- xxxx template align (for alignment) xxxxxxxxxxxxxxx -->

 <xsl:template name="align">
  <xsl:param name="min" />
  <xsl:param name="beginning" />

   <xsl:if test="$beginning &gt; $min">&#xA0;<xsl:call-template name="align">
       <xsl:with-param name="min" select="$min" />
       <xsl:with-param name="beginning" select="$beginning - 1" />

     </xsl:call-template>
   </xsl:if>
 </xsl:template>








<!-- xxxx  structure  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx -->

 <xsl:template match="structure">
  <h3> 
    <br />
    Structure: <xsl:value-of select="@name"/>
  </h3>
  <p>

    <xsl:if test="numGroups">
      Number of groups:
      <span>
        <xsl:value-of select="numGroups"/>
      </span>
      <br />
    </xsl:if>

    <xsl:if test="group">
      <xsl:for-each select="group">
        <xsl:number value="position()" format="1" />. group:
        <span>      
          <xsl:value-of select="."/>      
        </span>
        <br />
      </xsl:for-each>
    </xsl:if>

  </p>
 </xsl:template> 







<!-- xxxx  distance matrix  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx -->

 <xsl:template match="distanceMatrix">
   <h3>
     <br />
     Distance matrix: <xsl:value-of select="@name"/>
   </h3>  
   <p>

     <xsl:if test="matrixSize">
       Size of the distance matrix:
       <span>
         <xsl:value-of select="matrixSize"/>
       </span>
       <br />
     </xsl:if>

     Distance matrix:
   </p>

   <table>
     <xsl:call-template name="cr2br">
       <xsl:with-param name="text" select="matrix" />
       <xsl:with-param name="label" select="matrixLabels" />
     </xsl:call-template>

     <th class="without" />
     <xsl:call-template name="cr2Komma2">
       <xsl:with-param name="text" select="matrixLabels" />
     </xsl:call-template>

   </table>
   <br />
 </xsl:template>    



<!-- xxxx  cr2br: seperate string after "line break"  xxxxxxxxxxxxxx -->

 <xsl:template name="cr2br">
  <xsl:param name="text" />
  <xsl:param name="label" />

  <xsl:choose>
   <xsl:when test="contains($text, '&#xA;')"> 
    <tr>

      <xsl:choose>
        <xsl:when test="contains($label, ',')">
          <th class="without">
            <xsl:value-of select="substring-before($label, ',')" /> 
          </th>
        </xsl:when>
        <xsl:otherwise>
          <th class="without">
            <xsl:value-of select="$label" />
          </th>
        </xsl:otherwise>
      </xsl:choose>


      <xsl:call-template name="cr2Komma">
        <xsl:with-param name="text" select="substring-before($text, '&#xA;')" />
      </xsl:call-template>
    
    </tr>

    <xsl:call-template name="cr2br">
      <xsl:with-param name="text" select="substring-after($text, '&#xA;')" />
      <xsl:with-param name="label" select="substring-after($label, ',')" />
    </xsl:call-template>
   
   </xsl:when>
   <xsl:otherwise>
   
     <tr>

       <th class="without">
         <xsl:value-of select="$label" />
       </th>

       <xsl:call-template name="cr2Komma">
         <xsl:with-param name="text" select="$text" />
       </xsl:call-template>
     </tr>  

   </xsl:otherwise>
  </xsl:choose>
 </xsl:template>



<!-- xxxx  cr2br: seperate string after "comma"  xxxxxxxxxxxxxx -->

 <xsl:template name="cr2Komma">
  <xsl:param name="text" />

  <xsl:choose>
   <xsl:when test="contains($text, ',')">
    
     <td class="without">  
       <xsl:value-of select="substring-before($text, ',')" />
     </td>    

     <xsl:call-template name="cr2Komma">
       <xsl:with-param name="text" select="substring-after($text, ',')" />
     </xsl:call-template>
   </xsl:when>
   <xsl:otherwise>   
     
     <td class="without">
       <xsl:value-of select="$text" />
     </td>       

   </xsl:otherwise>
  </xsl:choose>
 </xsl:template>



<!-- xxxx  cr2br: seperate string after "comma" (for header)  xxxxxxxxxxxxxx -->

 <xsl:template name="cr2Komma2">
  <xsl:param name="text" />
  <xsl:choose>
   <xsl:when test="contains($text, ',')">


     <th class="without">
       <xsl:value-of select="substring-before($text, ',')" />
     </th>    

     <xsl:call-template name="cr2Komma2">
       <xsl:with-param name="text" select="substring-after($text, ',')" />
     </xsl:call-template>
   </xsl:when>
   <xsl:otherwise>   
     
     <th class="without">
       <xsl:value-of select="$text" />
     </th>       

   </xsl:otherwise>
  </xsl:choose>
 </xsl:template>



</xsl:stylesheet>