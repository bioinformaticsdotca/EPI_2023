-- callouts.lua
-- Handles static, regular/important/subtle dropdowns, and plain dropdowns.

local function get_default_icon(type)
  if type == "red" then return "fas fa-triangle-exclamation"
  elseif type == "blue" then return "fas fa-circle-info"
  elseif type == "green" then return "fas fa-circle-check"
  elseif type == "yellow" then return "fas fa-bell"
  elseif type == "purple" then return "fas fa-star"
  elseif type == "orange" then return "fas fa-fire"
  elseif type == "gray" then return "fas fa-comment-dots"
  else return "fas fa-info-circle"
  end
end

local function create_icon_element(icon_attr, type)
  local show_icon = true
  if icon_attr == "false" then show_icon = false end

  if show_icon then
    local icon_fa_class = (icon_attr and icon_attr ~= "true") and icon_attr or get_default_icon(type)
    -- Return the raw HTML for the icon
    return '<i class="' .. icon_fa_class .. '" aria-hidden="true"></i>'
  end
  return ''
end

local function create_title_element_str(title_text)
    if title_text and title_text ~= "" then
        local parsed_title_blocks = pandoc.read(title_text, 'markdown').blocks
        if parsed_title_blocks and #parsed_title_blocks > 0 then
            -- Use Pandoc to convert the title markdown to an HTML string
            return pandoc.write(pandoc.Pandoc(parsed_title_blocks), 'html'):gsub("^<p>", ""):gsub("</p>\n?$", "")
        else
            return pandoc.utils.stringify(title_text)
        end
    end
    return ''
end


function Div(el)
  if el.classes:includes("callout") then
    local type = el.attributes.type
    local style = el.attributes.style or "regular"
    local title_text = el.attributes.title
    local icon_attr = el.attributes.icon
    local center_title = el.attributes.center_title
    local is_collapsible = el.attributes.collapsible == "true"

    if not type or type == "" then
      io.stderr:write("Callout Warning: 'type' attribute is missing or empty for a .callout div.\n")
      return el
    end

    if is_collapsible then
      -- For collapsible callouts, we build the entire HTML as a single string to ensure correct nesting.
      local details_classes = {"callout-dropdown"}
      local base_type_class = "callout-" .. type

      if style == "plain" then
        table.insert(details_classes, "callout-plain")
        table.insert(details_classes, base_type_class)
      elseif style == "subtle" then
        table.insert(details_classes, "callout-subtle")
        table.insert(details_classes, base_type_class .. "-subtle")
      else -- regular or important
        table.insert(details_classes, base_type_class .. (style == "important" and "-important" or ""))
      end

      local summary_classes = {"callout-dropdown-summary"}
      local icon_html = create_icon_element(icon_attr, type)
      if icon_html ~= '' and style == "plain" then
        table.insert(summary_classes, "summary-has-icon")
      end

      local summary_title_text_val = title_text or "Details"
      local title_html = create_title_element_str(summary_title_text_val)

      local caret_html = '<span class="callout-dropdown-caret"><i class="fas fa-chevron-right" aria-hidden="true"></i></span>'

      -- Combine all parts of the summary, adding a non-breaking space after the icon.
      local summary_content_html = icon_html .. '&nbsp;&nbsp;&nbsp;<span class="callout-dropdown-title">' .. title_html .. '</span> ' .. caret_html
      local summary_html = '<summary class="'.. table.concat(summary_classes, " ") ..'">' .. summary_content_html .. '</summary>'

      -- Convert the div's content to HTML
      local body_content_html = pandoc.write(pandoc.Pandoc(el.content), 'html')
      local body_html = '<div class="callout-dropdown-content">' .. body_content_html .. '</div>'

      -- Assemble the final, complete <details> block
      local final_html = '<details class="' .. table.concat(details_classes, " ") .. '">' .. summary_html .. body_html .. '</details>'

      return pandoc.RawBlock('html', final_html)

    else -- Static callout (your original logic for this was sound)
      local specific_style_class = "callout-" .. type
      if style == "important" then specific_style_class = specific_style_class .. "-important"
      elseif style == "subtle" then specific_style_class = specific_style_class .. "-subtle"
      end

      -- create_title_element and create_icon_element need to return pandoc elements for the static case
      -- Re-defining local helper functions for the static case to return pandoc objects
      local function create_icon_element_static(icon_attr, type)
        local show_icon = true
        if icon_attr == "false" then show_icon = false end
        if show_icon then
          local icon_fa_class = (icon_attr and icon_attr ~= "true") and icon_attr or get_default_icon(type)
          local i_tag = pandoc.RawInline('html', '<i class="' .. icon_fa_class .. '" aria-hidden="true"></i>')
          return pandoc.Span({i_tag}, {class = "callout-icon"})
        end
        return nil
      end

      local function create_title_element_static(title_text, center_title_attr)
        if title_text and title_text ~= "" then
          local title_inline_content
          local parsed_title_blocks = pandoc.read(title_text, 'markdown').blocks
          if parsed_title_blocks and #parsed_title_blocks > 0 and parsed_title_blocks[1].t == "Para" then
              title_inline_content = parsed_title_blocks[1].content
          else
              title_inline_content = {pandoc.Str(title_text)}
          end
          local title_block_content = pandoc.Para(title_inline_content)
          local title_classes_list = {"callout-title-static"}
          if center_title_attr == "true" then table.insert(title_classes_list, "center") end
          return pandoc.Div({title_block_content}, {class=table.concat(title_classes_list, " ")})
        end
        return nil
      end

      local icon_element = create_icon_element_static(icon_attr, type)
      local title_element = create_title_element_static(title_text, center_title)

      local content_div_inner_blocks = {}
      if title_element then table.insert(content_div_inner_blocks, title_element) end
      for _, block in ipairs(el.content) do table.insert(content_div_inner_blocks, block) end
      local content_div_element = pandoc.Div(content_div_inner_blocks, {class = "callout-content"})

      local main_div_children = {}
      if icon_element then table.insert(main_div_children, icon_element) end
      table.insert(main_div_children, content_div_element)

      return pandoc.Div(main_div_children, {class = specific_style_class})
    end
  end
  return nil
end

-- A filter must return a list of filter definitions.
return { {Div = Div} }
