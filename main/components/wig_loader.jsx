/** @jsx html */
/// <reference no-default-lib="true"/>
/// <reference lib="dom" />
/// <reference lib="dom.asynciterable" />
/// <reference lib="deno.ns" />
import { html } from "../elemental.js"
import { Column, Row, Input, Code, EasyFilePicker, askForFiles } from "../elements.jsx"
import { data } from "../data.jsx"

// 
// file picker
// 
const CustomizedPicker = ({children, onChange})=><Row 
    style={`
        --default-width: 16.4rem;
        width: var(--default-width);
        max-width: var(--default-width);
        overflow-x: visible;
    `}
    >
        <EasyFilePicker
            style={`
                width: var(--default-width);
                max-width: var(--default-width);
                display: block;
                overflow-x: auto;
                background: #939393;
                position: relative;
                box-shadow: 0 4px 5px 0 rgba(0,0,0,0), 0 1px 10px 0 rgba(0,0,0,0), 0 2px 4px -1px rgba(0,0,0,0);
            `}
            hoverStyle={`
                min-width: max-content;
                overflow-x: visible;
                box-shadow: 0 4px 5px 0 rgba(0,0,0,0.14), 0 1px 10px 0 rgba(0,0,0,0.12), 0 2px 4px -1px rgba(0,0,0,0.3);
            `}
            onChange={onChange}
            margin-bottom="0.5rem"
            >
                <Row horizontalAlignment="center" width="100%" >
                    {children}
                </Row>
        </EasyFilePicker>
</Row>


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
            <CustomizedPicker onChange={(file)=>checkData({comWigFile: file})} >
                [Click to add Combined Wig]
            </CustomizedPicker>
            <CustomizedPicker onChange={(file)=>checkData({metadataFile: file})} >
                [Click to add Metadata]
            </CustomizedPicker>
    </Column>
}

// 
// column of all comwig+metadata pickers
// 
export const WigLoader = ({ children, style, }) => {
    let wigLoader
    // recursive picker (as soon as one is loaded add another)
    const addPicker = ()=>{
        wigLoader.appendChild(
            <CombinedWigElement 
                onLoaded={addPicker}
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
            min-width="21rem"
            max-width="21rem"
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
                        onLoaded={addPicker}
                        />
                </Column>}
        </Column>
    </Column>
}

