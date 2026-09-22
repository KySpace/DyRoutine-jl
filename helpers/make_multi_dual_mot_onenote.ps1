[CmdletBinding()]
param(
    [string]$DataRoot = 'C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data',
    [string]$OutputPath = '',
    [string]$PageTitle = 'Dual-isotope MOT table',
    [ValidateRange(20.0, 2000.0)]
    [double]$ImageDisplayWidth = 240.0,
    [switch]$ShowPage
)

# OneNote's desktop COM server is reliable through Windows PowerShell/.NET
# Framework, but its OpenHierarchy call fails through PowerShell 7 on this host.
if ($PSVersionTable.PSEdition -eq 'Core') {
    $windowsPowerShell = Join-Path $env:SystemRoot 'System32\WindowsPowerShell\v1.0\powershell.exe'
    if (-not (Test-Path -LiteralPath $windowsPowerShell -PathType Leaf)) {
        throw "Windows PowerShell is required for OneNote COM automation: $windowsPowerShell"
    }

    $relayArguments = @(
        '-NoProfile',
        '-ExecutionPolicy', 'Bypass',
        '-File', $PSCommandPath,
        '-DataRoot', $DataRoot,
        '-PageTitle', $PageTitle,
        '-ImageDisplayWidth', $ImageDisplayWidth.ToString([Globalization.CultureInfo]::InvariantCulture)
    )
    if (-not [string]::IsNullOrWhiteSpace($OutputPath)) {
        $relayArguments += @('-OutputPath', $OutputPath)
    }
    if ($ShowPage) {
        $relayArguments += '-ShowPage'
    }

    & $windowsPowerShell @relayArguments
    exit $LASTEXITCODE
}

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$script:OneNoteNamespace = 'http://schemas.microsoft.com/office/onenote/2013/onenote'
$script:ParentNames = @(
    'MOT loading 421',
    'MOT loading 626',
    'MOT lifetime',
    'CMOT lifetime',
    'ODT BField'
)
$script:PairNames = @(
    '162-164',
    '160-162',
    '161-162',
    '161-164',
    '163-164',
    '162-163',
    '161-163'
)

function Add-OneNoteElement {
    param(
        [Parameter(Mandatory)] [System.Xml.XmlDocument]$Document,
        [Parameter(Mandatory)] [System.Xml.XmlNode]$Parent,
        [Parameter(Mandatory)] [string]$Name,
        [hashtable]$Attributes = @{}
    )

    $node = $Document.CreateElement('one', $Name, $script:OneNoteNamespace)
    foreach ($attributeName in $Attributes.Keys) {
        $node.SetAttribute($attributeName, [string]$Attributes[$attributeName])
    }
    [void]$Parent.AppendChild($node)
    return $node
}

function Add-OneNoteText {
    param(
        [Parameter(Mandatory)] [System.Xml.XmlDocument]$Document,
        [Parameter(Mandatory)] [System.Xml.XmlNode]$Parent,
        [AllowEmptyString()] [string]$Text,
        [string]$Style = 'font-family:Calibri;font-size:10.0pt'
    )

    $oe = Add-OneNoteElement -Document $Document -Parent $Parent -Name 'OE' -Attributes @{ style = $Style }
    $textNode = Add-OneNoteElement -Document $Document -Parent $oe -Name 'T'
    [void]$textNode.AppendChild($Document.CreateCDataSection($Text))
    return $oe
}

function Add-OneNoteCell {
    param(
        [Parameter(Mandatory)] [System.Xml.XmlDocument]$Document,
        [Parameter(Mandatory)] [System.Xml.XmlNode]$Row,
        [string]$ShadingColor = ''
    )

    $attributes = @{}
    if ($ShadingColor) {
        $attributes.shadingColor = $ShadingColor
    }
    $cell = Add-OneNoteElement -Document $Document -Parent $Row -Name 'Cell' -Attributes $attributes
    return Add-OneNoteElement -Document $Document -Parent $cell -Name 'OEChildren'
}

function Get-PngDimensions {
    param([Parameter(Mandatory)] [byte[]]$Bytes, [Parameter(Mandatory)] [string]$Path)

    $signature = [byte[]](137, 80, 78, 71, 13, 10, 26, 10)
    if ($Bytes.Length -lt 24) {
        throw "PNG file is too short: $Path"
    }
    for ($index = 0; $index -lt $signature.Length; $index++) {
        if ($Bytes[$index] -ne $signature[$index]) {
            throw "Invalid PNG signature: $Path"
        }
    }

    $width = [uint32][System.Net.IPAddress]::NetworkToHostOrder(
        [BitConverter]::ToInt32($Bytes, 16)
    )
    $height = [uint32][System.Net.IPAddress]::NetworkToHostOrder(
        [BitConverter]::ToInt32($Bytes, 20)
    )
    if ($width -eq 0 -or $height -eq 0) {
        throw "PNG dimensions must be positive: $Path"
    }
    return [pscustomobject]@{ Width = [double]$width; Height = [double]$height }
}

function Get-PngEntries {
    param([Parameter(Mandatory)] [string]$Root)

    $entriesByCell = @{}
    $filenamePattern = '^\[([^\]]+)\]\.\[([^\]]+)\]\.\[([^\]]+)\]\.png$'

    foreach ($parentName in $script:ParentNames) {
        foreach ($pairName in $script:PairNames) {
            $cellKey = "$parentName`n$pairName"
            $entriesByCell[$cellKey] = [System.Collections.Generic.List[object]]::new()
            $pairPath = Join-Path (Join-Path $Root $parentName) $pairName
            if (-not (Test-Path -LiteralPath $pairPath -PathType Container)) {
                continue
            }

            foreach ($file in Get-ChildItem -LiteralPath $pairPath -File -Filter '*.png' | Sort-Object Name) {
                if ($file.Name -notmatch $filenamePattern) {
                    Write-Warning "Ignoring PNG with an unexpected name: $($file.FullName)"
                    continue
                }
                $bytes = [System.IO.File]::ReadAllBytes($file.FullName)
                $size = Get-PngDimensions -Bytes $bytes -Path $file.FullName
                $entriesByCell[$cellKey].Add([pscustomobject]@{
                    Path = $file.FullName
                    TestTag = $Matches[1]
                    StyleTag = $Matches[2]
                    SubvariantTag = $Matches[3]
                    Alt = "$($Matches[1]) / $($Matches[2]) / $($Matches[3])"
                    Width = $size.Width
                    Height = $size.Height
                    Base64 = [Convert]::ToBase64String($bytes)
                })
            }
        }
    }
    return $entriesByCell
}

function Get-ComparisonEntries {
    param([Parameter(Mandatory)] [string]$Root)

    $folder = Join-Path $Root 'Isotope pair comparison'
    $specs = @(
        @{ Key = '626-numbers'; Group = 'loading'; File = '[MOT.loading.pairs].[DCS-SCS].[nums].png' },
        @{ Key = '421-numbers'; Group = 'loading'; File = '[MOT.loading.pairs].[DDM-DIS].[nums].png' },
        @{ Key = '626-ratio'; Group = 'loading'; File = '[MOT.loading.pairs].[DCS-SCS].[ratio].png' },
        @{ Key = '421-ratio'; Group = 'loading'; File = '[MOT.loading.pairs].[DDM-DIS].[ratio].png' },
        @{ Key = 'cmot-values'; Group = 'lifetime'; File = '[CMOT.decay.pairs].[kappa].[values].png' },
        @{ Key = 'mot-values'; Group = 'lifetime'; File = '[MOT.decay.pairs].[tau].[values].png' },
        @{ Key = 'cmot-ratio'; Group = 'lifetime'; File = '[CMOT.decay.pairs].[DIS-DDM].[ratio].png' },
        @{ Key = 'mot-ratio'; Group = 'lifetime'; File = '[MOT.decay.pairs].[DDM-DIS].[ratio].png' }
    )
    $entries = @{}
    foreach ($spec in $specs) {
        $path = Join-Path $folder $spec.File
        if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
            throw "Missing isotope-pair comparison PNG: $path"
        }
        $bytes = [System.IO.File]::ReadAllBytes($path)
        $size = Get-PngDimensions -Bytes $bytes -Path $path
        $entries[$spec.Key] = [pscustomobject]@{
            Alt = "Isotope pair comparison / $($spec.Group) / $($spec.Key)"
            Width = $size.Width
            Height = $size.Height
            Base64 = [Convert]::ToBase64String($bytes)
        }
    }
    return $entries
}

function Add-OneNoteImage {
    param(
        [Parameter(Mandatory)] [System.Xml.XmlDocument]$Document,
        [Parameter(Mandatory)] [System.Xml.XmlNode]$Parent,
        [Parameter(Mandatory)] [psobject]$Entry,
        [Parameter(Mandatory)] [double]$DisplayWidth
    )

    $oe = Add-OneNoteElement -Document $Document -Parent $Parent -Name 'OE'
    $image = Add-OneNoteElement -Document $Document -Parent $oe -Name 'Image' -Attributes @{
        format = 'png'
        alt = $Entry.Alt
    }
    $displayHeight = $DisplayWidth * $Entry.Height / $Entry.Width
    [void](Add-OneNoteElement -Document $Document -Parent $image -Name 'Size' -Attributes @{
        width = $DisplayWidth.ToString('0.###', [Globalization.CultureInfo]::InvariantCulture)
        height = $displayHeight.ToString('0.###', [Globalization.CultureInfo]::InvariantCulture)
        isSetByUser = 'true'
    })
    $data = Add-OneNoteElement -Document $Document -Parent $image -Name 'Data'
    $data.InnerText = $Entry.Base64
}

function Add-ImageTable {
    param(
        [Parameter(Mandatory)] [System.Xml.XmlDocument]$Document,
        [Parameter(Mandatory)] [System.Xml.XmlNode]$CellChildren,
        [Parameter(Mandatory)] [System.Collections.IEnumerable]$Entries,
        [Parameter(Mandatory)] [double]$DisplayWidth,
        [Parameter(Mandatory)] [string]$CellLabel
    )

    $entryArray = @($Entries)
    if ($entryArray.Count -eq 0) {
        [void](Add-OneNoteText -Document $Document -Parent $CellChildren -Text '')
        return 0
    }

    $styles = @($entryArray.StyleTag | Sort-Object -Unique)
    $subvariants = @($entryArray.SubvariantTag | Sort-Object -Unique)
    if ($CellLabel -match '^(MOT loading|MOT lifetime|CMOT lifetime)') {
        $balanceOrder = @('t-balanced', 'n-balanced')
        $subvariants = @($balanceOrder | Where-Object { $_ -in $subvariants }) +
            @($subvariants | Where-Object { $_ -notin $balanceOrder })
    }
    $entryLookup = @{}
    foreach ($entry in $entryArray) {
        $key = "$($entry.StyleTag)`n$($entry.SubvariantTag)"
        if ($entryLookup.ContainsKey($key)) {
            throw "Duplicate style/subvariant position '$($entry.StyleTag)'/'$($entry.SubvariantTag)' in $CellLabel"
        }
        $entryLookup[$key] = $entry
    }

    $tableOe = Add-OneNoteElement -Document $Document -Parent $CellChildren -Name 'OE'
    $table = Add-OneNoteElement -Document $Document -Parent $tableOe -Name 'Table' -Attributes @{
        bordersVisible = 'true'
        hasHeaderRow = 'true'
    }

    $headerRow = Add-OneNoteElement -Document $Document -Parent $table -Name 'Row'
    $cornerChildren = Add-OneNoteCell -Document $Document -Row $headerRow -ShadingColor '#F2F2F2'
    [void](Add-OneNoteText -Document $Document -Parent $cornerChildren -Text '')
    foreach ($style in $styles) {
        $headerChildren = Add-OneNoteCell -Document $Document -Row $headerRow -ShadingColor '#F2F2F2'
        [void](Add-OneNoteText -Document $Document -Parent $headerChildren -Text $style -Style 'font-family:Calibri;font-size:9.0pt;font-weight:bold')
    }

    $imageCount = 0
    foreach ($subvariant in $subvariants) {
        $row = Add-OneNoteElement -Document $Document -Parent $table -Name 'Row'
        $labelChildren = Add-OneNoteCell -Document $Document -Row $row -ShadingColor '#F2F2F2'
        [void](Add-OneNoteText -Document $Document -Parent $labelChildren -Text $subvariant -Style 'font-family:Calibri;font-size:9.0pt;font-weight:bold')

        foreach ($style in $styles) {
            $imageChildren = Add-OneNoteCell -Document $Document -Row $row
            $key = "$style`n$subvariant"
            if (-not $entryLookup.ContainsKey($key)) {
                [void](Add-OneNoteText -Document $Document -Parent $imageChildren -Text '')
                continue
            }

            $entry = $entryLookup[$key]
            Add-OneNoteImage -Document $Document -Parent $imageChildren `
                -Entry $entry -DisplayWidth $DisplayWidth
            $imageCount++
        }
    }
    return $imageCount
}

function Add-ComparisonTable {
    param(
        [Parameter(Mandatory)] [System.Xml.XmlDocument]$Document,
        [Parameter(Mandatory)] [System.Xml.XmlNode]$OutlineChildren,
        [Parameter(Mandatory)] [hashtable]$Entries,
        [Parameter(Mandatory)] [double]$DisplayWidth
    )

    [void](Add-OneNoteText -Document $Document -Parent $OutlineChildren `
        -Text 'Isotope pair comparisons' `
        -Style 'font-family:Calibri;font-size:14.0pt;font-weight:bold')
    $groups = @(
        @{
            Caption = 'Loading numbers and ratios (626, 421)'
            Rows = @(
                @('626-numbers', '421-numbers'),
                @('626-ratio', '421-ratio')
            )
        },
        @{
            Caption = 'Lifetime parameters and ratios (1/κ, τ)'
            Rows = @(
                @('cmot-values', 'mot-values'),
                @('cmot-ratio', 'mot-ratio')
            )
        }
    )
    foreach ($group in $groups) {
        [void](Add-OneNoteText -Document $Document -Parent $OutlineChildren `
            -Text $group.Caption `
            -Style 'font-family:Calibri;font-size:10.0pt;font-weight:bold')
        $tableOe = Add-OneNoteElement -Document $Document -Parent $OutlineChildren -Name 'OE'
        $table = Add-OneNoteElement -Document $Document -Parent $tableOe -Name 'Table' -Attributes @{
            bordersVisible = 'true'
            hasHeaderRow = 'false'
        }
        foreach ($rowKeys in $group.Rows) {
            $row = Add-OneNoteElement -Document $Document -Parent $table -Name 'Row'
            foreach ($key in $rowKeys) {
                $cellChildren = Add-OneNoteCell -Document $Document -Row $row
                Add-OneNoteImage -Document $Document -Parent $cellChildren `
                    -Entry $Entries[$key] -DisplayWidth $DisplayWidth
            }
        }
    }
    return $Entries.Count
}

function New-OneNotePageXml {
    param(
        [Parameter(Mandatory)] [string]$PageId,
        [Parameter(Mandatory)] [string]$Title,
        [Parameter(Mandatory)] [hashtable]$EntriesByCell,
        [Parameter(Mandatory)] [hashtable]$ComparisonEntries,
        [Parameter(Mandatory)] [double]$DisplayWidth
    )

    $document = [System.Xml.XmlDocument]::new()
    $page = $document.CreateElement('one', 'Page', $script:OneNoteNamespace)
    $page.SetAttribute('ID', $PageId)
    $page.SetAttribute('name', $Title)
    [void]$document.AppendChild($page)

    $titleNode = Add-OneNoteElement -Document $document -Parent $page -Name 'Title'
    $titleOe = Add-OneNoteElement -Document $document -Parent $titleNode -Name 'OE' -Attributes @{
        style = 'font-family:Calibri;font-size:20.0pt'
    }
    $titleText = Add-OneNoteElement -Document $document -Parent $titleOe -Name 'T'
    [void]$titleText.AppendChild($document.CreateCDataSection($Title))

    $outline = Add-OneNoteElement -Document $document -Parent $page -Name 'Outline'
    [void](Add-OneNoteElement -Document $document -Parent $outline -Name 'Position' -Attributes @{ x = '36.0'; y = '90.0'; z = '0' })
    $outlineChildren = Add-OneNoteElement -Document $document -Parent $outline -Name 'OEChildren'
    $outerOe = Add-OneNoteElement -Document $document -Parent $outlineChildren -Name 'OE'
    $outerTable = Add-OneNoteElement -Document $document -Parent $outerOe -Name 'Table' -Attributes @{
        bordersVisible = 'true'
        hasHeaderRow = 'true'
    }

    $headerRow = Add-OneNoteElement -Document $document -Parent $outerTable -Name 'Row'
    $cornerChildren = Add-OneNoteCell -Document $document -Row $headerRow -ShadingColor '#D9EAF7'
    [void](Add-OneNoteText -Document $document -Parent $cornerChildren -Text 'Isotope pair' -Style 'font-family:Calibri;font-size:10.0pt;font-weight:bold')
    foreach ($parentName in $script:ParentNames) {
        $headerChildren = Add-OneNoteCell -Document $document -Row $headerRow -ShadingColor '#D9EAF7'
        [void](Add-OneNoteText -Document $document -Parent $headerChildren -Text $parentName -Style 'font-family:Calibri;font-size:10.0pt;font-weight:bold')
    }

    $imageCount = 0
    foreach ($pairName in $script:PairNames) {
        $row = Add-OneNoteElement -Document $document -Parent $outerTable -Name 'Row'
        $pairChildren = Add-OneNoteCell -Document $document -Row $row -ShadingColor '#D9EAF7'
        [void](Add-OneNoteText -Document $document -Parent $pairChildren -Text $pairName -Style 'font-family:Calibri;font-size:10.0pt;font-weight:bold')

        foreach ($parentName in $script:ParentNames) {
            $cellChildren = Add-OneNoteCell -Document $document -Row $row
            $cellKey = "$parentName`n$pairName"
            $imageCount += Add-ImageTable -Document $document -CellChildren $cellChildren `
                -Entries $EntriesByCell[$cellKey] -DisplayWidth $DisplayWidth `
                -CellLabel "$parentName/$pairName"
        }
    }

    $imageCount += Add-ComparisonTable -Document $document `
        -OutlineChildren $outlineChildren -Entries $ComparisonEntries `
        -DisplayWidth $DisplayWidth

    return [pscustomobject]@{ Document = $document; ImageCount = $imageCount }
}

if (-not (Test-Path -LiteralPath $DataRoot -PathType Container)) {
    throw "Data root does not exist: $DataRoot"
}
if ([string]::IsNullOrWhiteSpace($OutputPath)) {
    $OutputPath = Join-Path $DataRoot 'multi_dual_mot_table.one'
}
$OutputPath = [System.IO.Path]::GetFullPath($OutputPath)
$outputDirectory = Split-Path -Parent $OutputPath
if (-not (Test-Path -LiteralPath $outputDirectory -PathType Container)) {
    throw "Output directory does not exist: $outputDirectory"
}

$entriesByCell = Get-PngEntries -Root $DataRoot
$comparisonEntries = Get-ComparisonEntries -Root $DataRoot
$sourceImageCount = @($entriesByCell.Values | ForEach-Object { $_ }).Count + $comparisonEntries.Count
if ($sourceImageCount -eq 0) {
    throw "No matching PNG files were found under the configured parent/pair folders in $DataRoot"
}

$oneNote = $null
try {
    try {
        $oneNote = New-Object -ComObject OneNote.Application
    }
    catch {
        throw "Could not start the OneNote COM application: $($_.Exception.Message)"
    }

    $sectionId = ''
    try {
        $oneNote.OpenHierarchy($OutputPath, '', [ref]$sectionId, 3)
    }
    catch {
        throw "OneNote OpenHierarchy failed for '$OutputPath': $($_.Exception.Message) (HRESULT 0x$('{0:X8}' -f $_.Exception.HResult))"
    }

    $pageId = ''
    try {
        $oneNote.CreateNewPage($sectionId, [ref]$pageId, 0)
    }
    catch {
        throw "OneNote CreateNewPage failed for section '$sectionId': $($_.Exception.Message) (HRESULT 0x$('{0:X8}' -f $_.Exception.HResult))"
    }
    $page = New-OneNotePageXml -PageId $pageId -Title $PageTitle `
        -EntriesByCell $entriesByCell -ComparisonEntries $comparisonEntries `
        -DisplayWidth $ImageDisplayWidth
    if ($page.ImageCount -ne $sourceImageCount) {
        throw "Built $($page.ImageCount) image nodes from $sourceImageCount source images"
    }

    try {
        $oneNote.UpdatePageContent($page.Document.OuterXml, [datetime]::MinValue, 2, $true)
    }
    catch {
        throw "OneNote UpdatePageContent failed for page '$pageId': $($_.Exception.Message) (HRESULT 0x$('{0:X8}' -f $_.Exception.HResult))"
    }

    $verifiedXml = ''
    try {
        $oneNote.GetPageContent($pageId, [ref]$verifiedXml, 0, 2)
    }
    catch {
        throw "OneNote GetPageContent failed for page '$pageId': $($_.Exception.Message) (HRESULT 0x$('{0:X8}' -f $_.Exception.HResult))"
    }
    $verifiedDocument = [xml]$verifiedXml
    $namespaceManager = [System.Xml.XmlNamespaceManager]::new($verifiedDocument.NameTable)
    $namespaceManager.AddNamespace('one', $script:OneNoteNamespace)
    $verifiedImageCount = $verifiedDocument.SelectNodes('//one:Image', $namespaceManager).Count
    $verifiedTableCount = $verifiedDocument.SelectNodes('//one:Table', $namespaceManager).Count
    if ($verifiedImageCount -ne $sourceImageCount) {
        throw "OneNote returned $verifiedImageCount images after $sourceImageCount were submitted"
    }
    $verifiedComparisonCount = $verifiedDocument.SelectNodes(
        '//one:Image[starts-with(@alt, "Isotope pair comparison / ")]',
        $namespaceManager
    ).Count
    if ($verifiedComparisonCount -ne 8) {
        throw "OneNote returned $verifiedComparisonCount isotope-pair comparison images after 8 were submitted"
    }
    if ($verifiedTableCount -lt 3) {
        throw "OneNote returned only $verifiedTableCount table node(s)"
    }

    if ($ShowPage) {
        $oneNote.NavigateTo($pageId, '')
    }

    Write-Host "Created OneNote page '$PageTitle'."
    Write-Host "Section: $OutputPath"
    Write-Host "Embedded PNGs: $verifiedImageCount"
    Write-Host "Native tables: $verifiedTableCount"
    Write-Host "Page ID: $pageId"
}
finally {
    if ($null -ne $oneNote) {
        [void][Runtime.InteropServices.Marshal]::FinalReleaseComObject($oneNote)
    }
}
