function parse({csvString, seperator=",", useHeaders=null, firstRowIsHeaders=false, skipEmptyLines=true}) {
    let comments = []
    let rows = []
    let headers = []
    let isFirstDataRow = true

    lineLoop: for (const each of csvString.splitlines("\n")) {
        let eachLine = each.replace(/\r|\n/g, "")

        // 
        // comments
        // 
        if (eachLine.startsWith("#")) {
            comments.push(eachLine)
            continue
        }

        // 
        // empty lines
        // 
        if (skipEmptyLines && eachLine.trim().length == 0) {
            continue
        }

        let cells = eachLine.split(seperator)
        let cellsWithTypes = []
        let skipTo = 0
        let index = -1
        cellLoop: for (let eachCell of cells) {
            index += 1
            
            if (index < skipTo) {
                continue cellLoop
            }

            stripped = eachCell.trim()
            if (stripped.length == 0) {
                cellsWithTypes.push(null)
            } else {
                let firstChar = stripped[0]
                if ( !(firstChar.match(/"|\[|\{/)) ) {
                    try {
                        cellsWithTypes.push(JSON.parse(eachCell))
                    } catch (error) {
                        cellsWithTypes.push(eachCell)
                    }
                } else {
                    let remainingEndIndicies = []
                    let reverseIndex = cells.length
                    for (let _ of cells) {
                        reverseIndex--
                        remainingEndIndicies.push(reverseIndex)
                    }
                    for (const eachRemainingEndIndex of remainingEndIndicies) {
                        try {
                            const theSlice = cells.slice(index, eachRemainingEndIndex)
                            const partialString = theSlice.join(seperator)
                            cellsWithTypes.push(
                                JSON.parse(partialString)
                            )
                            continue cellLoop
                        } catch (error) {
                            
                        }
                    }
                    cellsWithTypes.push(eachCell)
                }
            }
        }

        // 
        // Headers
        // 
        if (isFirstDataRow) {
            isFirstDataRow = false
            if (firstRowIsHeaders) {
                headers = []
                for (let each of cellsWithTypes) {
                    headers.push(`${each}`)
                }
                continue lineLoop
            }
        }

        rows.push(cellsWithTypes)
    }

    if (firstRowIsHeaders || useHeaders) {
        // FIXME make each row have named items
    }

    return { comments, headers, rows }
    
}