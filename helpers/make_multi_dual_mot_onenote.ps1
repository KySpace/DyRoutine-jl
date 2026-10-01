[CmdletBinding()]
param(
    [string]$DataRoot = 'C:\Users\ky\OneDrive\Dy\DualIstpMOT\Data',
    [string]$OutputPath = '',
    [string]$PageTitle = 'Dual-isotope MOT table',
    [ValidateRange(20.0, 2000.0)]
    [double]$ImageDisplayWidth = 240.0,
    [switch]$CmotModelComparison,
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
    if ($CmotModelComparison) {
        $relayArguments += '-CmotModelComparison'
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
    'MOT loading balance',
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

function Get-ConfigSourceRecords {
    param([Parameter(Mandatory)] [string]$ConfigPath)

    $records = [System.Collections.Generic.List[object]]::new()
    $tag = ''
    $biases = @()
    foreach ($line in Get-Content -LiteralPath $ConfigPath) {
        if ($line -match '^\s*-\s*tag\s*:\s*"?([^"#]+)"?\s*$') {
            $tag = $Matches[1].Trim()
            $biases = @()
            continue
        }
        if ($line -match '^\s*-\s*tbiasmot\s*:\s*\[([^\]]*)\]') {
            $biases = @($Matches[1] -split ',' | ForEach-Object {
                $value = 0.0
                if ([double]::TryParse($_.Trim(),
                        [Globalization.NumberStyles]::Float,
                        [Globalization.CultureInfo]::InvariantCulture,
                        [ref]$value)) {
                    $value
                }
            })
            continue
        }
        if ($line -notmatch '^\s*source\s*:') {
            continue
        }
        foreach ($match in [regex]::Matches($line, 'liferes\s+(\d{4})\s+run(\d+)\.mat')) {
            $records.Add([pscustomobject]@{
                Tag = $tag
                Biases = @($biases)
                Date = $match.Groups[1].Value
                Run = [int]$match.Groups[2].Value
            })
        }
    }
    return $records
}

function Format-SourceCaption {
    param([Parameter(Mandatory)] [System.Collections.IEnumerable]$Records)

    $dateOrder = [System.Collections.Generic.List[string]]::new()
    $runsByDate = @{}
    foreach ($record in $Records) {
        if (-not $runsByDate.ContainsKey($record.Date)) {
            $dateOrder.Add($record.Date)
            $runsByDate[$record.Date] = [System.Collections.Generic.List[int]]::new()
        }
        if ($record.Run -notin $runsByDate[$record.Date]) {
            $runsByDate[$record.Date].Add($record.Run)
        }
    }
    $groups = foreach ($date in $dateOrder) {
        $runs = @($runsByDate[$date] | Sort-Object | ForEach-Object { $_.ToString('00') })
        if ($runs.Count -gt 0) {
            "$date-run$($runs -join ', ')"
        }
    }
    return $groups -join '; '
}

function Select-SourceRecords {
    param(
        [Parameter(Mandatory)] [System.Collections.IEnumerable]$Records,
        [Parameter(Mandatory)] [string]$SubvariantTag
    )

    $recordArray = @($Records)
    $selected = if ($SubvariantTag -eq 't-balanced') {
        @($recordArray | Where-Object { 0.0 -in $_.Biases })
    }
    elseif ($SubvariantTag -eq 'n-balanced') {
        @($recordArray | Where-Object {
            @($_.Biases | Where-Object { [Math]::Abs($_) -ge 1e-12 }).Count -gt 0
        })
    }
    elseif ($SubvariantTag -eq '*') {
        $recordArray
    }
    else {
        @($recordArray | Where-Object { $_.Tag -eq $SubvariantTag })
    }
    $selected = @($selected)
    return $(if ($selected.Count -eq 0) { $recordArray } else { $selected })
}

function Get-SourceCaption {
    param(
        [Parameter(Mandatory)] [string]$ConfigPath,
        [Parameter(Mandatory)] [string]$StyleTag,
        [Parameter(Mandatory)] [string]$SubvariantTag
    )

    if ($StyleTag -match 'fit' -or $StyleTag -match '(^|\.)log($|\.)') {
        return ''
    }
    $records = @(Get-ConfigSourceRecords -ConfigPath $ConfigPath)
    $selected = @(Select-SourceRecords -Records $records -SubvariantTag $SubvariantTag)
    return Format-SourceCaption -Records $selected
}

function Get-ParentSourceCaption {
    param(
        [Parameter(Mandatory)] [string]$Root,
        [Parameter(Mandatory)] [string]$ParentName,
        [Parameter(Mandatory)] [string]$SubvariantTag
    )

    $selected = [System.Collections.Generic.List[object]]::new()
    foreach ($pairName in $script:PairNames) {
        $configPath = Join-Path (Join-Path (Join-Path $Root $ParentName) $pairName) 'config.yaml'
        if (-not (Test-Path -LiteralPath $configPath -PathType Leaf)) {
            continue
        }
        $records = @(Get-ConfigSourceRecords -ConfigPath $configPath)
        if ($SubvariantTag -eq 'n-balanced' -and $pairName -ne '162-164' -and
                -not @($records | Where-Object {
                    @($_.Biases | Where-Object { [Math]::Abs($_) -ge 1e-12 }).Count -gt 0
                }).Count) {
            continue
        }
        foreach ($record in @(Select-SourceRecords -Records $records -SubvariantTag $SubvariantTag)) {
            $selected.Add($record)
        }
    }
    return Format-SourceCaption -Records $selected
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

            $configPath = Join-Path $pairPath 'config.yaml'
            foreach ($file in Get-ChildItem -LiteralPath $pairPath -File -Filter '*.png' | Sort-Object Name) {
                if ($file.Name -notmatch $filenamePattern) {
                    Write-Warning "Ignoring PNG with an unexpected name: $($file.FullName)"
                    continue
                }
                $testTag = $Matches[1]
                $styleTag = $Matches[2]
                $subvariantTag = $Matches[3]
                if ($parentName -eq 'MOT loading 421' -and
                    ($styleTag -notin @('lin', 'log', 'ratio.DDM-DIS', 'ratio.DIS-DCS') -or
                     -not $subvariantTag.Contains('.'))) {
                    continue
                }
                $bytes = [System.IO.File]::ReadAllBytes($file.FullName)
                $size = Get-PngDimensions -Bytes $bytes -Path $file.FullName
                $entriesByCell[$cellKey].Add([pscustomobject]@{
                    Path = $file.FullName
                    TestTag = $testTag
                    StyleTag = $styleTag
                    SubvariantTag = $subvariantTag
                    Alt = "$testTag / $styleTag / $subvariantTag"
                    Width = $size.Width
                    Height = $size.Height
                    Base64 = [Convert]::ToBase64String($bytes)
                    SourceCaption = if (Test-Path -LiteralPath $configPath -PathType Leaf) {
                        Get-SourceCaption -ConfigPath $configPath `
                            -StyleTag $styleTag -SubvariantTag $subvariantTag
                    } else { '' }
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
        @{ Key = 'odt-cmot-numbers'; Group = 'ODT/CMOT'; Parent = 'ODT BField'; Variant = '*'; File = '[ODT.CMOT.comparison].[log].[numbers].png' },
        @{ Key = 'odt-cmot-ratios'; Group = 'ODT/CMOT'; Parent = 'ODT BField'; Variant = '*'; File = '[ODT.CMOT.comparison].[lin].[ratios].png' },
        @{ Key = 't-421-numbers'; Group = 't-balanced'; Parent = 'MOT loading 421'; Variant = 't-balanced'; File = '[MOT.loading.pairs].[DDM-DIS-DCS].[nums.t-balanced].png' },
        @{ Key = 't-626-numbers'; Group = 't-balanced'; Parent = 'MOT loading 626'; Variant = '*'; File = '[MOT.loading.pairs].[DCS-SCS].[nums.t-balanced].png' },
        @{ Key = 't-cmot-values'; Group = 't-balanced'; Parent = 'CMOT lifetime'; Variant = 't-balanced'; File = '[CMOT.decay.pairs].[kappa].[values.t-balanced].png' },
        @{ Key = 't-cmot-n0'; Group = 't-balanced'; Parent = 'CMOT lifetime'; Variant = 't-balanced'; File = '[CMOT.decay.pairs].[n0].[values.t-balanced].png' },
        @{ Key = 't-mot-values'; Group = 't-balanced'; Parent = 'MOT lifetime'; Variant = 't-balanced'; File = '[MOT.decay.pairs].[tau].[values.t-balanced].png' },
        @{ Key = 't-mot-n0'; Group = 't-balanced'; Parent = 'MOT lifetime'; Variant = 't-balanced'; File = '[MOT.decay.pairs].[n0].[values.t-balanced].png' },
        @{ Key = 'n-421-numbers'; Group = 'n-balanced'; Parent = 'MOT loading 421'; Variant = 'n-balanced'; File = '[MOT.loading.pairs].[DDM-DIS].[nums.n-balanced].png' },
        @{ Key = 'n-cmot-values'; Group = 'n-balanced'; Parent = 'CMOT lifetime'; Variant = 'n-balanced'; File = '[CMOT.decay.pairs].[kappa].[values.n-balanced].png' },
        @{ Key = 'n-cmot-n0'; Group = 'n-balanced'; Parent = 'CMOT lifetime'; Variant = 'n-balanced'; File = '[CMOT.decay.pairs].[n0].[values.n-balanced].png' },
        @{ Key = 'n-mot-values'; Group = 'n-balanced'; Parent = 'MOT lifetime'; Variant = 'n-balanced'; File = '[MOT.decay.pairs].[tau].[values.n-balanced].png' },
        @{ Key = 'n-mot-n0'; Group = 'n-balanced'; Parent = 'MOT lifetime'; Variant = 'n-balanced'; File = '[MOT.decay.pairs].[n0].[values.n-balanced].png' },
        @{ Key = 't-421-ratio'; Group = 't-balanced'; Parent = 'MOT loading 421'; Variant = 't-balanced'; File = '[MOT.loading.pairs].[DDM-DIS].[ratio.t-balanced].png' },
        @{ Key = 't-421-dis-dcs-ratio'; Group = 't-balanced'; Parent = 'MOT loading 421'; Variant = 't-balanced'; File = '[MOT.loading.pairs].[DIS-DCS].[ratio.t-balanced].png' },
        @{ Key = 't-626-ratio'; Group = 't-balanced'; Parent = 'MOT loading 626'; Variant = '*'; File = '[MOT.loading.pairs].[DCS-SCS].[ratio.t-balanced].png' },
        @{ Key = 't-cmot-ratio'; Group = 't-balanced'; Parent = 'CMOT lifetime'; Variant = 't-balanced'; File = '[CMOT.decay.pairs].[DIS-DDM].[ratio.t-balanced].png' },
        @{ Key = 't-mot-ratio'; Group = 't-balanced'; Parent = 'MOT lifetime'; Variant = 't-balanced'; File = '[MOT.decay.pairs].[DDM-DIS].[ratio.t-balanced].png' },
        @{ Key = 'n-421-ratio'; Group = 'n-balanced'; Parent = 'MOT loading 421'; Variant = 'n-balanced'; File = '[MOT.loading.pairs].[DDM-DIS].[ratio.n-balanced].png' },
        @{ Key = 'n-cmot-ratio'; Group = 'n-balanced'; Parent = 'CMOT lifetime'; Variant = 'n-balanced'; File = '[CMOT.decay.pairs].[DIS-DDM].[ratio.n-balanced].png' },
        @{ Key = 'n-mot-ratio'; Group = 'n-balanced'; Parent = 'MOT lifetime'; Variant = 'n-balanced'; File = '[MOT.decay.pairs].[DDM-DIS].[ratio.n-balanced].png' }
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
            SourceCaption = Get-ParentSourceCaption -Root $Root `
                -ParentName $spec.Parent -SubvariantTag $spec.Variant
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

function Add-ImageWithCaption {
    param(
        [Parameter(Mandatory)] [System.Xml.XmlDocument]$Document,
        [Parameter(Mandatory)] [System.Xml.XmlNode]$Parent,
        [Parameter(Mandatory)] [psobject]$Entry,
        [Parameter(Mandatory)] [double]$DisplayWidth
    )

    Add-OneNoteImage -Document $Document -Parent $Parent `
        -Entry $Entry -DisplayWidth $DisplayWidth
    if ($Entry.PSObject.Properties.Name -contains 'SourceCaption' -and
            -not [string]::IsNullOrWhiteSpace($Entry.SourceCaption)) {
        [void](Add-OneNoteText -Document $Document -Parent $Parent `
            -Text $Entry.SourceCaption `
            -Style 'font-family:Calibri;font-size:7.0pt;text-align:center')
    }
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

    $styles = if ($CellLabel -like 'MOT loading 421/*') {
        @('lin', 'log', 'ratio.DDM-DIS', 'ratio.DIS-DCS')
    } elseif ($CellLabel -like 'MOT lifetime/*' -or $CellLabel -like 'CMOT lifetime/*') {
        @(@('lin', 'log', 'fit.log', 'size') | Where-Object { $_ -in $entryArray.StyleTag })
    } else {
        @($entryArray.StyleTag | Sort-Object -Unique)
    }
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
            Add-ImageWithCaption -Document $Document -Parent $imageChildren `
                -Entry $entry -DisplayWidth $DisplayWidth
            $imageCount++
        }
    }
    return $imageCount
}

function Add-ComparisonRows {
    param(
        [Parameter(Mandatory)] [System.Xml.XmlDocument]$Document,
        [Parameter(Mandatory)] [System.Xml.XmlNode]$Table,
        [Parameter(Mandatory)] [hashtable]$Entries,
        [Parameter(Mandatory)] [double]$DisplayWidth
    )

    $rows = @(
        @{ Label = 't-balanced values'; Cells = @('t-421-numbers', 't-626-numbers', $null, 't-mot-values', 't-cmot-values', 'odt-cmot-numbers') },
        @{ Label = 'n-balanced values'; Cells = @('n-421-numbers', $null, $null, 'n-mot-values', 'n-cmot-values', $null) },
        @{ Label = 't-balanced ratios'; Cells = @('t-421-ratio', 't-626-ratio', $null, 't-mot-ratio', 't-cmot-ratio', 'odt-cmot-ratios') },
        @{ Label = 'n-balanced ratios'; Cells = @('n-421-ratio', $null, $null, 'n-mot-ratio', 'n-cmot-ratio', $null) }
    )
    foreach ($rowSpec in $rows) {
        $row = Add-OneNoteElement -Document $Document -Parent $Table -Name 'Row'
        $labelChildren = Add-OneNoteCell -Document $Document -Row $row -ShadingColor '#D9EAF7'
        [void](Add-OneNoteText -Document $Document -Parent $labelChildren `
            -Text $rowSpec.Label -Style 'font-family:Calibri;font-size:10.0pt;font-weight:bold')
        foreach ($key in $rowSpec.Cells) {
            $cellChildren = Add-OneNoteCell -Document $Document -Row $row
            if ($null -eq $key) {
                [void](Add-OneNoteText -Document $Document -Parent $cellChildren -Text '')
            }
            else {
                Add-ImageWithCaption -Document $Document -Parent $cellChildren `
                    -Entry $Entries[$key] -DisplayWidth $DisplayWidth
                if ($key -eq 't-421-ratio') {
                    Add-ImageWithCaption -Document $Document -Parent $cellChildren `
                        -Entry $Entries['t-421-dis-dcs-ratio'] -DisplayWidth $DisplayWidth
                }
                $n0Key = switch ($key) {
                    't-cmot-values' { 't-cmot-n0' }
                    't-mot-values' { 't-mot-n0' }
                    'n-cmot-values' { 'n-cmot-n0' }
                    'n-mot-values' { 'n-mot-n0' }
                    default { $null }
                }
                if ($null -ne $n0Key) {
                    Add-ImageWithCaption -Document $Document -Parent $cellChildren `
                        -Entry $Entries[$n0Key] -DisplayWidth $DisplayWidth
                }
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

    $imageCount += Add-ComparisonRows -Document $document `
        -Table $outerTable -Entries $ComparisonEntries `
        -DisplayWidth $DisplayWidth

    return [pscustomobject]@{ Document = $document; ImageCount = $imageCount }
}

function New-CmotModelPageXmlAligned {
    param([Parameter(Mandatory)] [string]$PageId)

    $document = [System.Xml.XmlDocument]::new()
    $page = $document.CreateElement('one', 'Page', $script:OneNoteNamespace)
    $page.SetAttribute('ID', $PageId)
    $page.SetAttribute('name', 'CMOT decay model comparison')
    [void]$document.AppendChild($page)
    $titleNode = Add-OneNoteElement -Document $document -Parent $page -Name 'Title'
    $titleOe = Add-OneNoteElement -Document $document -Parent $titleNode -Name 'OE' -Attributes @{ style = 'font-family:Calibri;font-size:20.0pt' }
    $titleText = Add-OneNoteElement -Document $document -Parent $titleOe -Name 'T'
    [void]$titleText.AppendChild($document.CreateCDataSection('CMOT decay model comparison'))
    $outline = Add-OneNoteElement -Document $document -Parent $page -Name 'Outline'
    [void](Add-OneNoteElement -Document $document -Parent $outline -Name 'Position' -Attributes @{ x = '36.0'; y = '90.0'; z = '0' })
    $children = Add-OneNoteElement -Document $document -Parent $outline -Name 'OEChildren'
    $outerOe = Add-OneNoteElement -Document $document -Parent $children -Name 'OE'
    $table = Add-OneNoteElement -Document $document -Parent $outerOe -Name 'Table' -Attributes @{ bordersVisible = 'true'; hasHeaderRow = 'true' }
    $groupHeader = Add-OneNoteElement -Document $document -Parent $table -Name 'Row'
    foreach ($label in @('t-balanced', '', '', '', 'n-balanced', '', '', '')) {
        $cell = Add-OneNoteCell -Document $document -Row $groupHeader -ShadingColor '#B7D7EA'
        [void](Add-OneNoteText -Document $document -Parent $cell -Text $label -Style 'font-family:Calibri;font-size:9.0pt;font-weight:bold')
    }
    $columnHeader = Add-OneNoteElement -Document $document -Parent $table -Name 'Row'
    foreach ($label in @('Isotope pair', 'CMOT data', ':full', ':kappa', 'Isotope pair', 'CMOT data', ':full', ':kappa')) {
        $cell = Add-OneNoteCell -Document $document -Row $columnHeader -ShadingColor '#D9EAF7'
        [void](Add-OneNoteText -Document $document -Parent $cell -Text $label -Style 'font-family:Calibri;font-size:8.0pt;font-weight:bold')
    }
    $imageCount = 0
    foreach ($pairName in $script:PairNames) {
        $row = Add-OneNoteElement -Document $document -Parent $table -Name 'Row'
        foreach ($variant in @('t-balanced', 'n-balanced')) {
            $nameCell = Add-OneNoteCell -Document $document -Row $row -ShadingColor '#D9EAF7'
            [void](Add-OneNoteText -Document $document -Parent $nameCell -Text $pairName -Style 'font-family:Calibri;font-size:8.0pt;font-weight:bold')
            $pairFolder = Join-Path (Join-Path $DataRoot 'CMOT lifetime') $pairName
            $dataCell = Add-OneNoteCell -Document $document -Row $row
            $dataFiles = if (Test-Path -LiteralPath $pairFolder -PathType Container) {
                @(Get-ChildItem -LiteralPath $pairFolder -File -Filter '*.png' | Where-Object { $_.Name -match "^\[CMOT\.lifetime\]\.\[fit\.log\]\.\[$variant\]\.png$" } | Sort-Object Name)
            } else { @() }
            foreach ($file in $dataFiles) {
                $bytes = [System.IO.File]::ReadAllBytes($file.FullName)
                $size = Get-PngDimensions -Bytes $bytes -Path $file.FullName
                $entry = [pscustomobject]@{ Alt = "CMOT $variant data / $pairName"; Width = $size.Width; Height = $size.Height; Base64 = [Convert]::ToBase64String($bytes) }
                Add-OneNoteImage -Document $document -Parent $dataCell -Entry $entry -DisplayWidth 100.0
                $imageCount++
            }
            if (@($dataFiles).Count -eq 0) { [void](Add-OneNoteText -Document $document -Parent $dataCell -Text '') }
            foreach ($mode in @('full', 'kappa')) {
                $cell = Add-OneNoteCell -Document $document -Row $row
                $modeFolder = Join-Path $pairFolder 'alternatives'
                $modelFiles = if (Test-Path -LiteralPath $modeFolder -PathType Container) {
                    @(Get-ChildItem -LiteralPath $modeFolder -File -Filter '*.png' | Where-Object { $_.Name -match "^\[CMOT\.lifetime\]\.\[fit\.log\.$mode\]\.\[$variant\]\.png$" } | Sort-Object Name)
                } else { @() }
                foreach ($file in $modelFiles) {
                    $bytes = [System.IO.File]::ReadAllBytes($file.FullName)
                    $size = Get-PngDimensions -Bytes $bytes -Path $file.FullName
                    $entry = [pscustomobject]@{ Alt = "CMOT $variant $mode fit / $pairName"; Width = $size.Width; Height = $size.Height; Base64 = [Convert]::ToBase64String($bytes) }
                    Add-OneNoteImage -Document $document -Parent $cell -Entry $entry -DisplayWidth 100.0
                    $imageCount++
                }
                if (@($modelFiles).Count -eq 0) { [void](Add-OneNoteText -Document $document -Parent $cell -Text '') }
            }
        }
    }
    foreach ($parameter in @('n0', 'tau', 'kappa')) {
        $row = Add-OneNoteElement -Document $document -Parent $table -Name 'Row'
        $label = switch ($parameter) { 'n0' { 'N₀ comparison' } 'tau' { 'τ comparison' } default { 'κ comparison' } }
        foreach ($variant in @('t-balanced', 'n-balanced')) {
            $cell = Add-OneNoteCell -Document $document -Row $row -ShadingColor '#D9EAF7'
            [void](Add-OneNoteText -Document $document -Parent $cell -Text $label -Style 'font-family:Calibri;font-size:8.0pt;font-weight:bold')
            $cell = Add-OneNoteCell -Document $document -Row $row
            [void](Add-OneNoteText -Document $document -Parent $cell -Text '')
            foreach ($mode in @('full', 'kappa')) {
                $cell = Add-OneNoteCell -Document $document -Row $row
                if ($parameter -ne 'tau' -or $mode -eq 'full') {
                    $imagePath = Join-Path (Join-Path $DataRoot 'Isotope pair comparison') "[CMOT.models].[$parameter.$mode].[$variant].png"
                    if (-not (Test-Path -LiteralPath $imagePath -PathType Leaf)) { throw "Missing CMOT model comparison PNG: $imagePath" }
                    $bytes = [System.IO.File]::ReadAllBytes($imagePath)
                    $size = Get-PngDimensions -Bytes $bytes -Path $imagePath
                    $entry = [pscustomobject]@{ Alt = "CMOT $variant $parameter comparison / $mode"; Width = $size.Width; Height = $size.Height; Base64 = [Convert]::ToBase64String($bytes) }
                    Add-OneNoteImage -Document $document -Parent $cell -Entry $entry -DisplayWidth 100.0
                    $imageCount++
                } else { [void](Add-OneNoteText -Document $document -Parent $cell -Text '') }
            }
        }
    }
    $legendRow = Add-OneNoteElement -Document $document -Parent $table -Name 'Row'
    $legendCell = Add-OneNoteCell -Document $document -Row $legendRow
    $legendPath = Join-Path (Join-Path $DataRoot 'Isotope pair comparison') '[CMOT.models].[legend].[all].png'
    if (-not (Test-Path -LiteralPath $legendPath -PathType Leaf)) { throw "Missing CMOT model legend PNG: $legendPath" }
    $legendBytes = [System.IO.File]::ReadAllBytes($legendPath)
    $legendSize = Get-PngDimensions -Bytes $legendBytes -Path $legendPath
    $legendEntry = [pscustomobject]@{ Alt = 'CMOT model comparison marker legend'; Width = $legendSize.Width; Height = $legendSize.Height; Base64 = [Convert]::ToBase64String($legendBytes) }
    Add-OneNoteImage -Document $document -Parent $legendCell -Entry $legendEntry -DisplayWidth 120.0
    $imageCount++
    for ($index = 1; $index -lt 8; $index++) {
        $cell = Add-OneNoteCell -Document $document -Row $legendRow
        [void](Add-OneNoteText -Document $document -Parent $cell -Text '')
    }
    return [pscustomobject]@{ Document = $document; ImageCount = $imageCount }
}

function New-CmotModelPageXml {
    param([Parameter(Mandatory)] [string]$PageId)

    $document = [System.Xml.XmlDocument]::new()
    $page = $document.CreateElement('one', 'Page', $script:OneNoteNamespace)
    $page.SetAttribute('ID', $PageId)
    $page.SetAttribute('name', 'CMOT decay model comparison')
    [void]$document.AppendChild($page)
    $titleNode = Add-OneNoteElement -Document $document -Parent $page -Name 'Title'
    $titleOe = Add-OneNoteElement -Document $document -Parent $titleNode -Name 'OE' -Attributes @{
        style = 'font-family:Calibri;font-size:20.0pt'
    }
    $titleText = Add-OneNoteElement -Document $document -Parent $titleOe -Name 'T'
    [void]$titleText.AppendChild($document.CreateCDataSection('CMOT decay model comparison'))
    $outline = Add-OneNoteElement -Document $document -Parent $page -Name 'Outline'
    [void](Add-OneNoteElement -Document $document -Parent $outline -Name 'Position' -Attributes @{ x = '36.0'; y = '90.0'; z = '0' })
    $children = Add-OneNoteElement -Document $document -Parent $outline -Name 'OEChildren'
    $outerOe = Add-OneNoteElement -Document $document -Parent $children -Name 'OE'
    $table = Add-OneNoteElement -Document $document -Parent $outerOe -Name 'Table' -Attributes @{
        bordersVisible = 'true'
        hasHeaderRow = 'true'
    }
    $outerHeader = Add-OneNoteElement -Document $document -Parent $table -Name 'Row'
    foreach ($variant in @('t-balanced', 'n-balanced')) {
        $headerCell = Add-OneNoteCell -Document $document -Row $outerHeader -ShadingColor '#D9EAF7'
        [void](Add-OneNoteText -Document $document -Parent $headerCell -Text $variant -Style 'font-family:Calibri;font-size:11.0pt;font-weight:bold')
    }
    $outerRow = Add-OneNoteElement -Document $document -Parent $table -Name 'Row'
    $imageCount = 0
    foreach ($variant in @('t-balanced', 'n-balanced')) {
        $outerCell = Add-OneNoteCell -Document $document -Row $outerRow
        $innerOe = Add-OneNoteElement -Document $document -Parent $outerCell -Name 'OE'
        $innerTable = Add-OneNoteElement -Document $document -Parent $innerOe -Name 'Table' -Attributes @{
            bordersVisible = 'true'
            hasHeaderRow = 'true'
        }
        $innerHeader = Add-OneNoteElement -Document $document -Parent $innerTable -Name 'Row'
        foreach ($label in @('Isotope pair', 'CMOT data', ':full', ':kappa')) {
            $cell = Add-OneNoteCell -Document $document -Row $innerHeader -ShadingColor '#D9EAF7'
            [void](Add-OneNoteText -Document $document -Parent $cell -Text $label -Style 'font-family:Calibri;font-size:9.0pt;font-weight:bold')
        }
        foreach ($pairName in $script:PairNames) {
            $row = Add-OneNoteElement -Document $document -Parent $innerTable -Name 'Row'
            $nameCell = Add-OneNoteCell -Document $document -Row $row -ShadingColor '#D9EAF7'
            [void](Add-OneNoteText -Document $document -Parent $nameCell -Text $pairName -Style 'font-family:Calibri;font-size:8.0pt;font-weight:bold')
            $pairFolder = Join-Path (Join-Path $DataRoot 'CMOT lifetime') $pairName
            $dataCell = Add-OneNoteCell -Document $document -Row $row
            $dataFiles = if (Test-Path -LiteralPath $pairFolder -PathType Container) {
                @(Get-ChildItem -LiteralPath $pairFolder -File -Filter '*.png' |
                    Where-Object { $_.Name -match "^\[CMOT\.lifetime\]\.\[fit\.log\]\.\[$variant\]\.png$" } |
                    Sort-Object Name)
            } else { @() }
            foreach ($file in $dataFiles) {
                $bytes = [System.IO.File]::ReadAllBytes($file.FullName)
                $size = Get-PngDimensions -Bytes $bytes -Path $file.FullName
                $entry = [pscustomobject]@{ Alt = "CMOT $variant data / $pairName / $($file.BaseName)"; Width = $size.Width; Height = $size.Height; Base64 = [Convert]::ToBase64String($bytes) }
                Add-OneNoteImage -Document $document -Parent $dataCell -Entry $entry -DisplayWidth 125.0
                $imageCount++
            }
            if (@($dataFiles).Count -eq 0) { [void](Add-OneNoteText -Document $document -Parent $dataCell -Text '') }
            foreach ($mode in @('full', 'kappa')) {
                $cell = Add-OneNoteCell -Document $document -Row $row
                $modeFolder = Join-Path $pairFolder 'alternatives'
                $modelFiles = if (Test-Path -LiteralPath $modeFolder -PathType Container) {
                    @(Get-ChildItem -LiteralPath $modeFolder -File -Filter '*.png' |
                        Where-Object { $_.Name -match "^\[CMOT\.lifetime\]\.\[fit\.log\.$mode\]\.\[$variant\]\.png$" } |
                        Sort-Object Name)
                } else { @() }
                foreach ($file in $modelFiles) {
                    $bytes = [System.IO.File]::ReadAllBytes($file.FullName)
                    $size = Get-PngDimensions -Bytes $bytes -Path $file.FullName
                    $entry = [pscustomobject]@{ Alt = "CMOT $variant $mode fit / $pairName / $($file.BaseName)"; Width = $size.Width; Height = $size.Height; Base64 = [Convert]::ToBase64String($bytes) }
                    Add-OneNoteImage -Document $document -Parent $cell -Entry $entry -DisplayWidth 125.0
                    $imageCount++
                }
                if (@($modelFiles).Count -eq 0) { [void](Add-OneNoteText -Document $document -Parent $cell -Text '') }
            }
        }
        foreach ($parameter in @('n0', 'tau', 'kappa')) {
            $row = Add-OneNoteElement -Document $document -Parent $innerTable -Name 'Row'
            $labelCell = Add-OneNoteCell -Document $document -Row $row -ShadingColor '#D9EAF7'
            $parameterLabel = switch ($parameter) { 'n0' { 'N₀ comparison' } 'tau' { 'τ comparison' } default { 'κ comparison' } }
            [void](Add-OneNoteText -Document $document -Parent $labelCell -Text $parameterLabel -Style 'font-family:Calibri;font-size:8.0pt;font-weight:bold')
            $emptyCell = Add-OneNoteCell -Document $document -Row $row
            [void](Add-OneNoteText -Document $document -Parent $emptyCell -Text '')
            foreach ($mode in @('full', 'kappa')) {
                $cell = Add-OneNoteCell -Document $document -Row $row
                $available = $parameter -ne 'tau' -or $mode -eq 'full'
                if ($available) {
                    $imagePath = Join-Path (Join-Path $DataRoot 'Isotope pair comparison') "[CMOT.models].[$parameter.$mode].[$variant].png"
                    if (-not (Test-Path -LiteralPath $imagePath -PathType Leaf)) { throw "Missing CMOT model comparison PNG: $imagePath" }
                    $bytes = [System.IO.File]::ReadAllBytes($imagePath)
                    $size = Get-PngDimensions -Bytes $bytes -Path $imagePath
                    $entry = [pscustomobject]@{ Alt = "CMOT $variant $parameter comparison / $mode"; Width = $size.Width; Height = $size.Height; Base64 = [Convert]::ToBase64String($bytes) }
                    Add-OneNoteImage -Document $document -Parent $cell -Entry $entry -DisplayWidth 125.0
                    $imageCount++
                }
                else { [void](Add-OneNoteText -Document $document -Parent $cell -Text '') }
            }
        }
        $imageCount += 0
    }
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

if ($CmotModelComparison) {
    $oneNote = $null
    try {
        try { $oneNote = New-Object -ComObject OneNote.Application }
        catch { throw "Could not start the OneNote COM application: $($_.Exception.Message)" }
        $sectionId = ''
        $oneNote.OpenHierarchy($OutputPath, '', [ref]$sectionId, 3)
        $pageId = ''
        $oneNote.CreateNewPage($sectionId, [ref]$pageId, 0)
        $page = New-CmotModelPageXmlAligned -PageId $pageId
        $oneNote.UpdatePageContent($page.Document.OuterXml, [datetime]::MinValue, 2, $true)
        $verifiedXml = ''
        $oneNote.GetPageContent($pageId, [ref]$verifiedXml, 0, 2)
        $verifiedDocument = [xml]$verifiedXml
        $namespaceManager = [System.Xml.XmlNamespaceManager]::new($verifiedDocument.NameTable)
        $namespaceManager.AddNamespace('one', $script:OneNoteNamespace)
        $verifiedImages = $verifiedDocument.SelectNodes('//one:Image', $namespaceManager).Count
        if ($verifiedImages -ne $page.ImageCount) {
            throw "OneNote returned $verifiedImages model images after $($page.ImageCount) were submitted"
        }
        $verifiedTables = $verifiedDocument.SelectNodes('//one:Table', $namespaceManager)
        if ($verifiedTables.Count -ne 1) { throw "Expected one aligned CMOT comparison table; OneNote returned $($verifiedTables.Count) tables" }
        $modelTable = $verifiedTables[0]
        $columnCount = $modelTable.SelectNodes('./one:Row[1]/one:Cell', $namespaceManager).Count
        $rowCount = $modelTable.SelectNodes('./one:Row', $namespaceManager).Count
        if ($columnCount -ne 8 -or $rowCount -ne 13) {
            throw "CMOT comparison table has $columnCount columns and $rowCount rows; expected 8 and 13"
        }
        if ($ShowPage) { $oneNote.NavigateTo($pageId, '') }
        Write-Host "Created OneNote page 'CMOT decay model comparison'."
        Write-Host "Section: $OutputPath"
        Write-Host "Embedded PNGs: $verifiedImages"
        Write-Host "Aligned balance table: 8 columns"
        Write-Host "Page ID: $pageId"
    }
    finally {
        if ($null -ne $oneNote) { [void][Runtime.InteropServices.Marshal]::FinalReleaseComObject($oneNote) }
    }
    return
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
    if ($verifiedComparisonCount -ne $comparisonEntries.Count) {
        throw "OneNote returned $verifiedComparisonCount isotope-pair comparison images after $($comparisonEntries.Count) were submitted"
    }
    if ($verifiedTableCount -lt 3) {
        throw "OneNote returned only $verifiedTableCount table node(s)"
    }

    $outerTable = $verifiedDocument.SelectSingleNode('//one:Table', $namespaceManager)
    foreach ($pairName in @('161-163', '162-164')) {
        $pairIndex = [Array]::IndexOf($script:PairNames, $pairName) + 1
        $pairRow = $outerTable.SelectNodes('./one:Row', $namespaceManager)[$pairIndex]
        $loadingTable = $pairRow.SelectSingleNode('./one:Cell[2]//one:Table', $namespaceManager)
        if ($null -eq $loadingTable) {
            throw "OneNote read-back is missing the 421 table for $pairName"
        }
        $columnCount = $loadingTable.SelectNodes('./one:Row[1]/one:Cell', $namespaceManager).Count - 1
        $rowCount = $loadingTable.SelectNodes('./one:Row', $namespaceManager).Count - 1
        $expectedRows = @($EntriesByCell["MOT loading 421`n$pairName"] |
            ForEach-Object { $_.SubvariantTag } | Sort-Object -Unique).Count
        if ($columnCount -ne 4 -or $rowCount -ne $expectedRows) {
            throw "OneNote 421 table for $pairName has $columnCount figure columns and $rowCount tag/balance rows; expected 4 and $expectedRows"
        }
        Write-Host "Verified OneNote $pairName 421 table: $columnCount figure columns, $rowCount tag/balance rows."
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
