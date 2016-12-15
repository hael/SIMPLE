<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0"
                xmlns="http://www.w3.org/1999/xhtml"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                xmlns:cml="http://www.xml-cml.org/schema">

<xsl:output indent="yes"/>

<xsl:template match="/">
    <xsl:apply-templates select="/cml:dictionary|/cml:cml/cml:dictionary" />
</xsl:template>

<xsl:template match="//cml:dictionary">
<html>
<head>
    <title>CML Dictionary - <xsl:value-of select="@title" /></title>
    <meta http-equiv="Content-Type" content="application/xhtml+xml; charset=utf-8" />
    <link rel="stylesheet" href="http://www.xml-cml.org/common/style/cml-dict.css" type="text/css" />
</head>

<body>
<div class="head">
    <h1><xsl:value-of select="@title" /></h1>

    <h2>Namespace</h2>
    <p>The namespace of this dictionary is: <code><xsl:value-of select="@namespace" /></code></p>

<xsl:if test="@dictionaryPrefix">
    <h2>Default Prefix</h2>
    <p>The default prefix for this dictionary is: <code><xsl:value-of select="@dictionaryPrefix" /></code></p>
</xsl:if>

<xsl:if test="cml:description">
    <h2>Description</h2>
    <div><xsl:copy-of select="cml:description"/></div>
</xsl:if>

    <hr />
</div>

<h1>Table of Contents</h1>

<ul>
<xsl:for-each select="./cml:entry">
    <li><a href="#{@id}"><xsl:value-of select="@term" /></a></li>
</xsl:for-each>
</ul>

<hr />

<xsl:apply-templates select="./cml:entry"/>

</body>
</html>

</xsl:template>

<xsl:template match="//cml:entry">

    <div class="entry">

    <h1 title="{@id}"><a name="{@id}"><xsl:value-of select="@term"/> (ID: <xsl:value-of select="@id" />)</a></h1>
    <xsl:if test="cml:definition">
        <h2>Definition</h2>
        <div><xsl:copy-of select="cml:definition"/></div>
    </xsl:if>
    <xsl:if test="cml:description">
        <h2>Description</h2>
        <div><xsl:copy-of select="cml:description"/></div>
    </xsl:if>
    <xsl:if test="@dataType">
        <h2>Data Type</h2>
        <p><em><xsl:value-of select="@term"/></em> is of data type <code><xsl:value-of select="@dataType"/></code></p>
    </xsl:if>
    <xsl:if test="@unitType">
        <h2>Unit Type</h2>
        <p><em><xsl:value-of select="@term"/></em> has unit type <code><xsl:value-of select="@unitType"/></code></p>
    </xsl:if>
    <xsl:if test="@unit">
        <h2>Default units</h2>
        <p><em><xsl:value-of select="@term"/></em> has default units <code><xsl:value-of select="@unit"/></code></p>
    </xsl:if>

    <hr />

    </div>

</xsl:template>

</xsl:stylesheet>