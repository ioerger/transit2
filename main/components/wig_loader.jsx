/** @jsx html */
/// <reference no-default-lib="true"/>
/// <reference lib="dom" />
/// <reference lib="dom.asynciterable" />
/// <reference lib="deno.ns" />
import { html } from "../elemental.js"
import { Column, Row, Input, Code, askForFiles } from "../elements.jsx"
import { Event, trigger, everyTime, once } from "https://deno.land/x/good@0.5.15/events.js"
import { EasyFilePicker } from "./file_picker.jsx"
import { data, events } from "../data.jsx"


// 
// Comwig+Metadata picker
// 
const CombinedWigElement = ({onLoaded})=>{
    const data = {
        comWigFile: null,
        metadataFile: null,
    }
    // helper that knows when to call "onLoaded"
    let onLoadedWasCalled = false
    const checkData = (newData)=>{
        // integrate the new data
        Object.assign(data, newData)
        console.log(`wig data is:`,data)
        // check
        if (!onLoadedWasCalled && data.comWigFile && data.metadataFile) {
            onLoadedWasCalled = true
            onLoaded(data)
        }
    }
    return <Column
        style={`
            padding: 1rem;
            border: gray 1px solid;
            border-radius: 1rem;
            background: rgba(255,255,255,0.6);
            box-shadow: 0 4px 5px 0 rgba(0,0,0,0.14), 0 1px 10px 0 rgba(0,0,0,0.12), 0 2px 4px -1px rgba(0,0,0,0.3);
            margin-bottom: 1rem;
        `}
        >
            <EasyFilePicker defaultWidth="17rem" onChange={(file)=>checkData({comWigFile: file})} >
                [Click to add Combined Wig]
            </EasyFilePicker>
            <EasyFilePicker defaultWidth="17rem" onChange={(file)=>checkData({metadataFile: file})} >
                [Click to add Metadata]
            </EasyFilePicker>
    </Column>
}

// 
// column of all comwig+metadata pickers
// 
export const WigLoader = ({ children, style, }) => {
    let wigLoader
    // recursive picker (as soon as one is loaded add another)
    const onLoaded = (newWigGroup)=>{
        data.wigGroups.push(newWigGroup)
        console.log(`events.wigDataAdded`)
        trigger(events.wigDataAdded)
        // add another picker for additional comwig files
        wigLoader.appendChild(
            <CombinedWigElement 
                onLoaded={onLoaded}
                />
        )
    }
    
    return <Column name="WigLoader" height="100%">
        <span class="custom-header">
            Wig Files
        </span>
        <Column
            name="WigLoader-Column"
            padding="1.2rem 1rem"
            horizontalAlignment={`space-between`}
            min-width="21.5rem"
            max-width="21.5rem"
            background="rgba(50, 134, 153, 0.2)"
            border="lightgray solid 1px"
            border-radius="0.7rem"
            height="100%"
            overflow="visible"
            >
                {/* a Wrapper to make scrolling work and overflow: visible work */}
                {wigLoader = <Column
                    hoverStyle={`
                        min-width: 150vw;
                    `}
                    overflow="auto"
                    height="100%"
                    >
                    <CombinedWigElement 
                        onLoaded={onLoaded}
                        />
                </Column>}
        </Column>
    </Column>
}

