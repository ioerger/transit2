/** @jsx html */
/// <reference no-default-lib="true"/>
/// <reference lib="dom" />
/// <reference lib="dom.asynciterable" />
/// <reference lib="deno.ns" />
import { html } from "../elemental.js"
import { Event, trigger, everyTime, once } from "https://deno.land/x/good@0.5.15/events.js"
import { Column, Row, Input, Code, askForFiles } from "../elements.jsx"
import { EasyFilePicker } from "./file_picker.jsx"
import { data, events } from "../data.jsx"

const standardColumnNames = [ "disabled", "name", "condition", ]

const createGridRowElements = ({ rowData, columns }) => {
    console.debug(`rowData is:`,rowData)
    let rowElements = [
        <Input type="checkbox" checked={rowData.disabled} onChange={(event)=>{ rowData.disabled = event.target.checked }} />
    ]
    for (const eachColumn of columns) {
        // skip disabled column, thats done manually
        if (eachColumn == "disabled") {
            continue
        }
        const value = rowData[eachColumn]
        // if primitive, like string or number
        if (value != null && !(rowData[eachColumn] instanceof Object)) {
            rowElements.push(<div class="custom-grid-cell"> {rowData[eachColumn]} </div>)
        } else {
            rowElements.push(<div class="custom-grid-cell"> </div>)
        }
    }
    return rowElements
}

export const GridElement = ({columnNames, children})=> {
    return <Column
            padding="1.2rem 1rem"
            background="rgba(187, 225, 161, 0.28)"
            border="lightgray solid 1px"
            border-radius="0.7rem"
            flex-grow="1"
            style={`
                display: grid;
                grid-template-columns: ${columnNames.map(each=>"auto").join(" ")};
                width: 100%;
                grid-column-gap: 1px;
                grid-row-gap: 5px;
            `}
            >
                {children}
        </Column>
}

export const SampleFileTable = ({ children, style, }) => {
    let sampleFileElement, gridElement
    // when new data is added
    everyTime(events.wigDataAdded).then(()=>{
        console.log(`everyTime(events.wigDataAdded)`)
        const columns = new Set(standardColumnNames)
        // create all the columns first
        for (const eachSample of data.samples) {
            for (const [key, value] of Object.entries(eachSample)) {
                if (!(value instanceof Object)) {
                    columns.add(key)
                }
            }
        }

        const elements = []
        const headerElements = []
        for (const eachColumn of columns) {
            headerElements.push(
                <Row background="gray" color="white" padding="0.4rem 0.8rem">
                    {eachColumn}
                </Row>
            )
        }
        elements.push(headerElements)
        for (const eachSample of data.samples) {
            elements.push(createGridRowElements({ columns, rowData: eachSample }))
        }
        
        // create the new grid
        const newGrid = <GridElement columnNames={[...columns]}>
            {elements.flat()}
        </GridElement>
        // replace the old grid
        gridElement.remove()
        sampleFileElement.appendChild(newGrid)
        gridElement = newGrid
    })
    
    return sampleFileElement = <Column flex-grow="1" width="100%">
        <span class="custom-header">
            Sample Files
        </span>
        {gridElement = <GridElement columnNames={standardColumnNames} />}
    </Column>
}

