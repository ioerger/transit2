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

const MenuItem = ({children, text, onClick}) => {
    const childContainer = <div style={`position: absolute; right: 0; top: -2px; display: none; transform: translateX(100%); flex-direction: column; border-radius: 5px; z-index: 100;`}>
        {children}
    </div>
    let container
    return container = <Column
        class="custom-menu-container"
        style={`position: relative; padding: 0.7rem; background: whitesmoke; border-left: 2px solid gray; border: 2px solid gray;`}
        onMouseOver={()=>{
            childContainer.style.display = "flex"
            container.style.background = "white"
        }}
        onMouseOut={()=>{
            childContainer.style.display = "none"
            container.style.background = "whitesmoke"
        }}
        >
            <Row onClick={onClick} style={`width: max-content;`}>
                {text}
            </Row>
            {childContainer}
    </Column>
}

export const Menu = ({ children, style, }) => {
    return <Row style="background: whitesmoke;">
        <MenuItem text="File">
            <MenuItem text="Export">
                <MenuItem text="Selected Datasets"></MenuItem>
            </MenuItem>
            <MenuItem text="Convert">
                <MenuItem text="prot_table to PTT"> </MenuItem>
                <MenuItem text="prot_table to GFF3"> </MenuItem>
            </MenuItem>
        </MenuItem>
        <MenuItem text="View">
        </MenuItem>
        <MenuItem text="Analysis">
            <MenuItem text="Himar1 Methods">
                <MenuItem text="Gumble" />
                <MenuItem text="Resampling" />
            </MenuItem>
            <MenuItem text="Tn5 Methods">
                <MenuItem text="Gumble" />
                <MenuItem text="Resampling" />
            </MenuItem>
        </MenuItem>
        <MenuItem text="Help">
        </MenuItem>
    </Row>
}

